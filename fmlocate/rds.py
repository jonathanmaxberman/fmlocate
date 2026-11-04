"""Minimal RDS/RBDS decoder operating on FM multiplex (MPX) audio.

Copied from RadioShelf. Input is the demodulated MPX signal, for example as produced by
``rtl_fm -M fm -l 0 -A std -p 0 -s 171k -F 9`` (signed 16-bit mono, 171 kHz).
171 kHz / 24 = 7125 Hz = exactly 6 samples per RDS bit (1187.5 bit/s).
"""
from __future__ import annotations

from collections import Counter, defaultdict

import numpy as np
from scipy import signal

FS = 171_000
DECIM = 24
SPB = 6  # samples per bit after decimation

# Offset words for block A, B, C, C', D
OFFSETS = {0x0FC: "A", 0x198: "B", 0x168: "C", 0x350: "Cp", 0x1B4: "D"}
POLY = 0x5B9  # x^10 + x^8 + x^7 + x^5 + x^4 + x^3 + 1

PTY_RBDS = [
    "None", "News", "Information", "Sports", "Talk", "Rock", "Classic Rock",
    "Adult Hits", "Soft Rock", "Top 40", "Country", "Oldies", "Soft", "Nostalgia",
    "Jazz", "Classical", "Rhythm and Blues", "Soft R&B", "Foreign Language",
    "Religious Music", "Religious Talk", "Personality", "Public", "College",
    "Spanish Talk", "Spanish Music", "Hip Hop", "", "", "Weather", "Emergency Test",
    "Emergency",
]


def _check(info: int) -> int:
    reg = info << 10
    for i in range(25, 9, -1):
        if reg & (1 << i):
            reg ^= POLY << (i - 10)
    return reg & 0x3FF


_CHECK_TABLE = np.array([_check(i) for i in range(1 << 16)], dtype=np.uint16)


def _syndrome(word: int) -> int:
    return int(_CHECK_TABLE[word >> 10]) ^ (word & 0x3FF)


# Error patterns (bursts up to 3 bits) keyed by syndrome, for correcting blocks.
_BURSTS: dict[int, int] = {}
for _span, _fills in ((1, [0b1]), (2, [0b11]), (3, [0b101, 0b111])):
    for _pat in _fills:
        for _shift in range(0, 27 - _span):
            _BURSTS.setdefault(_syndrome(_pat << _shift), _pat << _shift)

BLOCK_OFFSETS = [(0x0FC,), (0x198,), (0x168, 0x350), (0x1B4,)]


def _baseband(mpx: np.ndarray) -> np.ndarray:
    """RDS complex baseband at 7125 Hz, carrier locked to 3x the stereo pilot when present."""
    n = np.arange(len(mpx))
    w19 = 2 * np.pi * 19_000 / FS
    lp_pilot = signal.butter(2, 50, fs=FS, output="sos")
    pilot = signal.sosfiltfilt(lp_pilot, mpx * np.exp(-1j * w19 * n))
    amp = np.abs(pilot)
    use_pilot = np.median(amp) > 4 * np.std(amp) / np.sqrt(len(amp) / 1000) and \
        np.median(amp) > 0.3 * np.sqrt(np.mean(mpx ** 2)) * 0.05
    sos = signal.butter(6, 2_400, fs=FS, output="sos")
    if use_pilot:
        ph = 3 * np.unwrap(np.angle(pilot))
        bb = signal.sosfiltfilt(sos, mpx * np.exp(-1j * (3 * w19 * n + ph)))[::DECIM]
        theta = np.angle(np.sum(bb * bb)) / 2  # constant pilot/RDS phase offset
        return bb * np.exp(-1j * theta)
    bb = signal.sosfiltfilt(sos, mpx * np.exp(-2j * np.pi * 57_000 / FS * n))[::DECIM]
    sq = bb * bb
    win = 2 * SPB * 64
    k = np.ones(win) / win
    sm = np.convolve(sq.real, k, "same") + 1j * np.convolve(sq.imag, k, "same")
    return bb * np.exp(-1j * np.unwrap(np.angle(sm)) / 2)


def _bits_from_mpx(mpx: np.ndarray) -> np.ndarray:
    """Recover differential-decoded RDS data bits from MPX samples."""
    bb = _baseband(mpx)
    if len(bb) < SPB * 200:
        return np.zeros(0, dtype=np.uint8)
    r = bb.real
    shape = np.r_[np.ones(SPB // 2), -np.ones(SPB // 2)]
    mf = np.convolve(r, shape[::-1], "full")[SPB - 1:]
    # Bit clock is locked to the pilot and the dongle has a TCXO, so a single
    # timing phase per ~2 s window is accurate; choose it by maximum energy.
    bitwin = 2048 * SPB
    positions = []
    for start in range(0, len(mf) - SPB, bitwin):
        seg = mf[start:start + bitwin]
        nb = len(seg) // SPB
        if nb == 0:
            break
        energy = [np.abs(seg[o:o + nb * SPB:SPB]).mean() for o in range(SPB)]
        o = int(np.argmax(energy))
        positions.extend(range(start + o, start + nb * SPB, SPB))
    pos = np.array([p for p in positions if p < len(mf)])
    sym = (mf[pos] > 0).astype(np.uint8)
    return (sym[1:] ^ sym[:-1]).astype(np.uint8)


def _words(bits: np.ndarray) -> np.ndarray:
    weights = (1 << np.arange(25, -1, -1)).astype(np.int64)
    win = np.lib.stride_tricks.sliding_window_view(bits.astype(np.int64), 26)
    return win @ weights


def _blocks(bits: np.ndarray):
    """(position, type, info) for every bit position holding an error-free block."""
    if len(bits) < 26:
        return []
    words = _words(bits)
    info = (words >> 10).astype(np.int64)
    offs = _CHECK_TABLE[info].astype(np.int64) ^ (words & 0x3FF)
    out = []
    for off, name in OFFSETS.items():
        for p in np.nonzero(offs == off)[0]:
            out.append((int(p), name, int(info[p])))
    out.sort()
    return out


def _correct(word: int, offsets: tuple) -> int | None:
    for off in offsets:
        syn = _syndrome(word ^ off)
        if syn == 0:
            return word >> 10
        pat = _BURSTS.get(syn)
        if pat is not None:
            return (word ^ pat) >> 10
    return None


def _groups(bits: np.ndarray):
    """Yield groups as lists of 4 block infos (None where uncorrectable)."""
    for _, g in groups_at(bits):
        yield g


def groups_at(bits: np.ndarray):
    """Yield (bit position, [4 block infos]) for each group; position / 1187.5 = seconds."""
    blocks = _blocks(bits)
    idx = {"A": 0, "B": 1, "C": 2, "Cp": 2, "D": 3}
    # Group start positions implied by error-free blocks, voted per phase.
    starts = sorted({p - 26 * idx[t] for p, t, _ in blocks if p - 26 * idx[t] >= 0})
    if not starts:
        return
    phase_votes = Counter(s % 104 for s in starts)
    good = {s for s in starts if phase_votes[s % 104] >= 2}
    if not good:
        return
    words = _words(bits)
    seq = sorted(good)
    full = set(seq)
    for a, b in zip(seq, seq[1:]):  # fill gaps between aligned groups
        if (b - a) % 104 == 0:
            full.update(range(a + 104, b, 104))
    for g in sorted(full):
        if g + 104 > len(bits):
            continue
        yield g, [_correct(int(words[g + 26 * i]), BLOCK_OFFSETS[i]) for i in range(4)]


def decode(mpx: np.ndarray) -> dict:
    """Decode RDS from MPX samples. Returns dict with pi, ps, rt, pty, stats."""
    mpx = np.asarray(mpx, dtype=np.float64)
    mpx = mpx - mpx.mean()
    bits = _bits_from_mpx(mpx)

    pi_votes: Counter = Counter()
    ps_votes = [Counter() for _ in range(8)]
    rt_votes: dict[int, Counter] = defaultdict(Counter)
    pty_votes: Counter = Counter()
    groups = 0
    for a, b, c, d in _groups(bits):
        if a is None or b is None:
            continue
        groups += 1
        pi_votes[a] += 1
        gtype = b >> 12
        version_b = (b >> 11) & 1
        pty_votes[(b >> 5) & 0x1F] += 1
        if gtype == 0 and d is not None:
            seg = b & 0x3
            for i, ch in enumerate((d >> 8, d & 0xFF)):
                if 32 <= ch < 127:
                    ps_votes[seg * 2 + i][chr(ch)] += 1
        elif gtype == 2:
            seg = b & 0xF
            if not version_b and c is not None and d is not None:
                for i, ch in enumerate((c >> 8, c & 0xFF, d >> 8, d & 0xFF)):
                    rt_votes[seg * 4 + i][ch] += 1
            elif version_b and d is not None:
                for i, ch in enumerate((d >> 8, d & 0xFF)):
                    rt_votes[seg * 2 + i][ch] += 1

    result = {
        "groups": groups,
        "bits": int(len(bits)),
        "pi": None, "callsign": None, "ps": None, "rt": None, "pty": None,
    }
    if not pi_votes:
        return result
    pi, n = pi_votes.most_common(1)[0]
    # Require a clear majority so random corrected noise can't produce an ID.
    if n < 3 or n < 0.5 * groups:
        return result
    result["pi"] = f"{pi:04X}"
    result["callsign"] = pi_to_callsign(pi)
    # Only trust a station name whose every character was received clearly;
    # scrolling/rotating PS text would otherwise produce a mix of messages.
    clear = all(v and v.most_common(1)[0][1] >= 2 and
                v.most_common(1)[0][1] >= 0.6 * sum(v.values()) for v in ps_votes)
    if clear:
        result["ps"] = "".join(v.most_common(1)[0][0] for v in ps_votes).strip() or None
    if rt_votes:
        chars = []
        for i in range(64):
            if i not in rt_votes:
                chars.append(" ")
                continue
            ch = rt_votes[i].most_common(1)[0][0]
            if ch == 0x0D:
                break
            chars.append(chr(ch) if 32 <= ch < 127 else " ")
        rt = " ".join("".join(chars).split())
        result["rt"] = rt or None
    if pty_votes:
        result["pty"] = PTY_RBDS[pty_votes.most_common(1)[0][0]]
    return result


def pi_to_callsign(pi: int) -> str | None:
    """Convert a North American RBDS PI code to a call sign (NRSC-4)."""
    if pi >> 12 == 0xA:
        pi = (((pi >> 8) & 0xF) << 12) | (pi & 0xFF)
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if 4096 <= pi <= 21671:
        prefix, n = "K", pi - 4096
    elif 21672 <= pi <= 39247:
        prefix, n = "W", pi - 21672
    else:
        return None
    return prefix + letters[n // 676] + letters[(n // 26) % 26] + letters[n % 26]


def decode_file(path: str) -> dict:
    return decode(np.fromfile(path, dtype="<i2"))


if __name__ == "__main__":
    import json
    import sys
    print(json.dumps(decode_file(sys.argv[1]), indent=2))
