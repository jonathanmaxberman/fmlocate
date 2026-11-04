"""Live measurement: sweep the FM band, decode RDS on the strongest channels, then sweep
the TV bands and the NOAA Weather Radio channels."""
from __future__ import annotations

import json
import os
import time
import urllib.request
from datetime import datetime, timedelta, timezone

import numpy as np
from scipy import signal

from .fcc import TV_CHANNELS, tv_channel_low_mhz
from .noaa import CHANNELS_MHZ as WX_CHANNELS
from . import rds
from .solve import ip_location

# RadioShelf (a local radio recorder) shares the dongle; if it's running, don't interrupt it.
RADIOSHELF_API = os.environ.get("FMLOCATE_RADIOSHELF_URL", "http://localhost:8765/api")
CENTERS = np.round(np.arange(87.9e6, 108.0e6, 200e3))
TV_BANDS = {"vlo": (53e6, 89e6), "vhi": (173e6, 217e6), "uhf": (469e6, 609e6)}
ATSC_PILOT_HZ = 309_441          # above the channel's lower edge
ATSC_PILOT_FRACTION = 10 ** (-11.3 / 10)  # pilot power / total ATSC 1.0 signal power
ATSC_BW_HZ = 5.381e6
PILOT_MIN_DB = 4.0               # 3-bin pilot excess counted as a detection
PILOT_MAX_DB = 14.0              # above this the pilot no longer tracks S/N (saturation, bin errors)
POWER_MIN_DB = 6.0               # channel power above floor counted as a detection without a pilot


class TunerBusy(RuntimeError):
    pass


def check_tuner_free(margin_min: float = 15) -> str:
    """Refuse to touch the dongle if RadioShelf is using it or will need it soon.
    Returns None when RadioShelf isn't running (the usual case)."""
    try:
        with urllib.request.urlopen(f"{RADIOSHELF_API}/status", timeout=3) as r:
            status = json.load(r)
    except (OSError, ValueError):
        return None
    if status.get("recording") or (status.get("tuner") or {}).get("owner"):
        raise TunerBusy(f"RadioShelf is using the tuner: {status.get('tuner')}")
    if (status.get("scan") or {}).get("running"):
        raise TunerBusy("RadioShelf is scanning")
    with urllib.request.urlopen(f"{RADIOSHELF_API}/upcoming", timeout=5) as r:
        items = json.load(r)
    if isinstance(items, dict):
        items = items.get("items", [])
    now = datetime.now(timezone.utc)
    soon = now + timedelta(minutes=margin_min)
    for it in items:
        start = datetime.fromisoformat(it["start"])
        end = datetime.fromisoformat(it["end"])
        if start < soon and end > now:
            raise TunerBusy(f"Recording '{it['title']}' due at {it['start']}")
    return "RadioShelf idle, no recording within %d min" % margin_min


def channel_levels(f: np.ndarray, p: np.ndarray, half_width: float = 60e3,
                   floor_window: float = 1.5e6, floor_pct: float = 20) -> list[dict]:
    """Mean dB within +-half_width of each channel, relative to a local noise floor.

    RadioShelf uses the band-wide median as the floor; the floor here rises ~8 dB from
    88 to 108 MHz, so it is the 20th percentile of the bins within +-1.5 MHz instead
    (mostly the gaps between stations).
    """
    out = []
    for c in CENTERS:
        m = np.abs(f - c) < half_width
        w = np.abs(f - c) < floor_window
        if m.any():
            floor = float(np.percentile(p[w], floor_pct))
            out.append({"freq": round(c / 1e6, 1), "level_db": round(float(p[m].mean() - floor), 2)})
    return out


def tv_band_of(ch: int) -> str:
    return "vlo" if ch <= 6 else "vhi" if ch <= 13 else "uhf"


def tv_channel_levels(bands: dict) -> list[dict]:
    """Per-channel S/N estimate (dB, over the 5.38 MHz signal) from averaged TV-band spectra.

    Two measures, because neither works across the whole range:
    - Pilot: ATSC 1.0 carries a CW pilot 11.3 dB below its total power. In a ~5 kHz bin it
      stands ~19 dB above the data, so its excess over the in-channel bins detects stations
      down to about -17 dB S/N without knowing the noise floor. It saturates once the data
      outweighs the noise (S/N above ~0 dB), and multipath can notch it.
    - Power: mean channel power over the band's floor (25th percentile of channel means).
      Used once it is clearly above the floor, and the only measure for ATSC 3.0 (no pilot).
    """
    out = []
    for name, b in bands.items():
        f = np.array(b["f_hz"], dtype=float)
        P = 10 ** (np.array(b["db"]) / 10)
        chans = [c for c in TV_CHANNELS if tv_band_of(c) == name]
        lo = {c: tv_channel_low_mhz(c) * 1e6 for c in chans}
        mean = {c: float(np.mean(P[(f > lo[c] + 0.5e6) & (f < lo[c] + 5.5e6)])) for c in chans}
        floor = float(np.percentile(list(mean.values()), 25))
        binw = float(np.median(np.diff(f)))
        k_pilot = ATSC_PILOT_FRACTION * ATSC_BW_HZ / binw
        # rtl_power's bin labels sit about one bin off; find the shift from clear pilots, then
        # look within +-1.5 bins of it.
        def pilot(c, centre, half):
            pc = lo[c] + ATSC_PILOT_HZ + centre
            w = np.flatnonzero(np.abs(f - pc) <= half)
            if len(w) == 0:
                return None, 0.0
            k = int(w[np.argmax(P[w])])
            ref = np.median(P[(f > lo[c] + 0.4e6) & (f < lo[c] + 0.9e6)])
            ex = (P[max(k - 1, 0):k + 2].sum() - 3 * ref) / ref
            return float(10 * np.log10(1 + max(ex, -0.99))), float(f[k] - lo[c] - ATSC_PILOT_HZ)
        first = {c: pilot(c, 3e3, 10e3) for c in chans}
        offs = [o for c, (x, o) in first.items() if x is not None and x > 10]
        shift = float(np.median(offs)) if offs else 3e3
        for c in chans:
            x, _ = pilot(c, shift, 1.5 * binw)
            x = x if x is not None else 0.0
            level = 10 * np.log10(mean[c] / floor)
            q = (10 ** (min(x, PILOT_MAX_DB) / 10) - 1) / k_pilot
            snr_pilot = 10 * np.log10(q / (1 - q)) if 0 < q < 0.9 else None
            snr_power = 10 * np.log10(max(10 ** (level / 10) - 1, 1e-3))
            skip = False
            if level >= POWER_MIN_DB and x < PILOT_MIN_DB and name != "uhf":
                # Local electrical noise is strong on VHF and looks like this; can't tell.
                how, snr, skip = "power without pilot (noise?)", None, True
            elif level >= POWER_MIN_DB:
                how, snr = ("power" if x >= PILOT_MIN_DB else "power, no pilot (ATSC 3.0?)"), snr_power
            elif x >= PILOT_MIN_DB and snr_pilot is not None:
                how, snr = "pilot", snr_pilot
            else:
                how, snr = None, None
            out.append({"channel": c, "freq": tv_channel_low_mhz(c) + 3.0, "band": name,
                        "level_db": round(float(level), 2), "pilot_db": round(x, 2), "how": how,
                        "snr_db": None if snr is None else round(float(snr), 2), "skip": skip,
                        "thr_db": round(tv_threshold_snr_db(binw), 2)})
    return out


def tv_threshold_snr_db(bin_hz: float = 5.3e3) -> float:
    """S/N at which the pilot just reaches PILOT_MIN_DB: the censoring threshold for TV."""
    q = (10 ** (PILOT_MIN_DB / 10) - 1) / (ATSC_PILOT_FRACTION * ATSC_BW_HZ / bin_hz)
    return float(10 * np.log10(q / (1 - q)))


def wx_channel_levels(spec: dict) -> list[dict]:
    """NOAA Weather Radio levels: mean dB within +-7 kHz minus the 20th percentile of the span."""
    f = np.array(spec["f_hz"], dtype=float)
    p = np.array(spec["db"])
    floor = 10 ** (np.percentile(p, 20) / 10)
    out = []
    for c in WX_CHANNELS:
        m = np.abs(f - c * 1e6) <= 7e3
        out.append({"freq": c, "level_db": round(float(10 * np.log10(np.mean(10 ** (p[m] / 10)) / floor)), 2)})
    return out


def levels_from_spectrum(meas: dict) -> list[dict]:
    s = meas["spectrum"]
    p = np.array(s["db"])
    f = np.array(s["f_hz"]) if "f_hz" in s else s["f0_hz"] + np.arange(len(p)) * s["df_hz"]
    return channel_levels(f, p)


def _retry(fn, ok, tries: int = 3, wait: float = 7, log=print):
    """Retry captures that come back empty (another program briefly holding the device, e.g.
    RadioShelf's rtl_test probe when its status is polled)."""
    for i in range(tries):
        out = fn()
        if ok(out):
            return out
        if i < tries - 1:
            log(f"(device busy, retrying in {wait:.0f} s)", end=" ")
            time.sleep(wait)
    return out


def pilot_snr_db(mpx: np.ndarray, fs: int = 171_000) -> float:
    f, p = signal.welch(mpx.astype(float), fs, nperseg=8192)
    pilot = p[np.abs(f - 19_000) < 30].max()
    noise = np.median(p[(f > 60_000) & (f < 75_000)])
    return float(10 * np.log10(pilot / (noise + 1e-12)))


def _averaged(sdr, lo, hi, step, seconds, sweeps, label, log):
    """Several sweeps averaged in power, on the first sweep's frequency grid."""
    spectra = []
    for i in range(sweeps):
        log(f"Sweeping {label} ({i + 1}/{sweeps}, {seconds} s)...")
        f, p = _retry(lambda: sdr.sweep(lo, hi, step, seconds), lambda r: len(r[0]) > 100, log=log)
        if not len(f):
            raise RuntimeError("the sweep returned no data (is the SDR attached and free?)")
        if spectra:
            p = np.interp(f0, f, p)
        else:
            f0 = f
        spectra.append(p)
    p = 10 * np.log10(np.mean([10 ** (s / 10) for s in spectra], axis=0))
    return {"f_hz": [int(round(x)) for x in f0], "db": [round(float(x), 2) for x in p]}


def measure_tv(sdr, seconds: int = 5, sweeps: int = 4, log=print) -> dict:
    bands = {name: _averaged(sdr, lo, hi, 10e3, seconds, sweeps,
                             f"TV {name} {lo / 1e6:.0f}-{hi / 1e6:.0f} MHz", log)
             for name, (lo, hi) in TV_BANDS.items()}
    return {"seconds": seconds, "sweeps": sweeps, "bands": bands}


def measure_wx(sdr, seconds: int = 10, sweeps: int = 2, log=print) -> dict:
    spec = _averaged(sdr, 161.9e6, 163.0e6, 1e3, seconds, sweeps, "NOAA 162.4-162.55 MHz", log)
    return {"seconds": seconds, "sweeps": sweeps, "spectrum": spec}


def measure(sdr_spec: str = "auto", gain: float | None = None, ppm: int = 0, sweep_seconds: int = 10,
            sweeps: int = 6, probe_min_db: float = 3.0, probe_max: int = 14, probe_seconds: float = 8,
            extra: tuple = ("tv", "wx"), log=print) -> dict:
    from .sdr import open_backend
    status = check_tuner_free()
    if status:
        log(f"Tuner: {status}")
        time.sleep(7)  # our /api/status call may have started RadioShelf's rtl_test probe
    sdr = open_backend(sdr_spec, gain, ppm, log)
    try:
        spectra = []
        for i in range(sweeps):
            log(f"Sweeping 87.5-108.3 MHz ({i + 1}/{sweeps}, {sweep_seconds} s)...")
            f, p = _retry(lambda: sdr.sweep(87.5e6, 108.3e6, 10e3, sweep_seconds),
                          lambda r: len(r[0]) > 1000, log=log)
            if not len(f):
                raise RuntimeError("the sweep returned no data (is the SDR attached and free?)")
            if spectra:
                p = np.interp(f0, f, p)
            else:
                f0 = f
            spectra.append(p)
        # Average in the power domain (indoor fading moves levels by ~10 dB over minutes).
        f = f0
        p = 10 * np.log10(np.mean([10 ** (s / 10) for s in spectra], axis=0))
        per_sweep = [channel_levels(f, s) for s in spectra]
        levels = channel_levels(f, p)
        lv = {c["freq"]: c["level_db"] for c in levels}
        cands = [c for c in levels if c["level_db"] >= probe_min_db
                 and max(lv.get(round(c["freq"] + d, 1), -99) for d in (-0.2, 0.2)) - c["level_db"] <= 10]
        cands = sorted(cands, key=lambda c: -c["level_db"])[:probe_max]
        probes = []
        for i, c in enumerate(cands):
            log(f"RDS {c['freq']:.1f} MHz ({c['level_db']:.1f} dB) [{i + 1}/{len(cands)}]...", end=" ")
            mpx = _retry(lambda: sdr.capture_mpx(c["freq"] * 1e6, probe_seconds),
                         lambda m: len(m) >= 171_000, log=log)
            if len(mpx) < 171_000:
                log("no audio")
                continue
            info = rds.decode(mpx)
            pilot = pilot_snr_db(mpx)
            log(f"pilot {pilot:.0f} dB, PI {info['pi']} {info['callsign'] or ''} {info['ps'] or ''}")
            probes.append({"freq": c["freq"], "pi": info["pi"], "callsign": info["callsign"],
                           "ps": info["ps"], "pilot_db": round(pilot, 1)})
        tv = measure_tv(sdr, log=log) if "tv" in extra else None
        wx = measure_wx(sdr, log=log) if "wx" in extra else None
    finally:
        sdr.close()
    return {
        "time_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "receiver": {"sdr": sdr.describe(), "gain_db": sdr.gain, "ppm": ppm, "sweep_seconds": sweep_seconds, "sweeps": sweeps,
                     "crop": 0.3,
                     "level_definition": "mean dB over +-60 kHz of channel minus 20th percentile "
                                         "of bins within +-1.5 MHz"},
        "per_sweep_levels": [[c["level_db"] for c in lv] for lv in per_sweep],
        # Saved so a replay can fall back to the same coarse prior without the network.
        "ip_location": ip_location(),
        "channels": levels,
        "rds": probes,
        "spectrum": {"f_hz": [int(round(x)) for x in f],
                     "db": [round(float(x), 2) for x in p]},
        **({"tv": tv} if tv else {}),
        **({"wx": wx} if wx else {}),
    }
