"""Receiver backends. Each one can sweep a frequency range into a power spectrum and
capture demodulated FM multiplex (for RDS).

- `rtl`:   the rtl-sdr command-line tools (rtl_power, rtl_fm). RTL-SDR dongles only.
- `soapy`: SoapySDR's Python bindings, which cover most SDRs (RTL-SDR, HackRF, Airspy,
           SDRplay, LimeSDR, PlutoSDR, USRP, ...) through their Soapy modules.

`open_backend("auto")` uses the rtl tools if they find a dongle, else the first
SoapySDR device.
"""
from __future__ import annotations

import os
import select
import shutil
import subprocess
import tempfile
import time
from fractions import Fraction

import numpy as np
from scipy import signal

EXTRA_PATHS = ["/opt/homebrew/bin", "/usr/local/bin"]
MPX_FS = 171_000  # what rds.py expects
FM_IF_FS = 240_000


class NoDevice(RuntimeError):
    pass


def binary(name: str) -> str | None:
    path = os.pathsep.join([os.environ.get("PATH", "")] + EXTRA_PATHS)
    return shutil.which(name, path=path)


def fm_demod_mpx(iq: np.ndarray, fs: float) -> np.ndarray:
    """Complex baseband centred on the station -> FM multiplex at 171 kHz (float)."""
    fr = Fraction(FM_IF_FS / fs).limit_denominator(2000)
    x = signal.resample_poly(iq, fr.numerator, fr.denominator)
    mpx = np.angle(x[1:] * np.conj(x[:-1]))
    return signal.resample_poly(mpx, 57, 80).astype(np.float32) * 10_000


class RtlBackend:
    name = "rtl"

    def __init__(self, gain: float | None = None, ppm: int = 0, index: int = 0):
        if not binary("rtl_power") or not binary("rtl_fm"):
            raise NoDevice("rtl_power/rtl_fm not found (macOS: brew install librtlsdr; "
                           "Debian/Ubuntu: apt install rtl-sdr)")
        self.gain = 40.2 if gain is None else gain
        self.ppm, self.index = ppm, index

    def describe(self) -> str:
        return f"rtl-sdr tools, device {self.index}, gain {self.gain} dB"

    @staticmethod
    def probe(index: int = 0) -> str | None:
        """Name of the attached dongle, or None."""
        exe = binary("rtl_test")
        if not exe:
            return None
        try:
            p = subprocess.run([exe, "-d", str(index), "-t"], capture_output=True, text=True, timeout=8)
            out = p.stdout + p.stderr
        except subprocess.TimeoutExpired as e:
            out = (e.stdout or b"").decode(errors="replace") + (e.stderr or b"").decode(errors="replace")
        if "No supported devices" in out:
            return None
        for line in out.splitlines():
            s = line.strip()
            if s[:2].strip().rstrip(":").isdigit() and ":" in s[:4]:
                return s.split(":", 1)[1].strip()
        return None

    def sweep(self, f_lo, f_hi, bin_hz, seconds, crop=0.3):
        """One rtl_power pass. -c drops each hop's filter roll-off, which otherwise leaves
        5-8 dB steps in the floor and under-reads channels near hop boundaries."""
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
            path = tmp.name
        try:
            cmd = [binary("rtl_power"), "-d", str(self.index), "-f", f"{int(f_lo)}:{int(f_hi)}:{int(bin_hz)}",
                   "-g", str(self.gain), "-p", str(self.ppm), "-c", str(crop), "-i", str(int(seconds)),
                   "-1", path]
            try:
                subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                               timeout=seconds * 12 + 60)
            except subprocess.TimeoutExpired:
                pass
            f, p = [], []
            with open(path) as fh:
                for line in fh:
                    r = line.strip().split(",")
                    if len(r) < 7:
                        continue
                    lo, st = float(r[2]), float(r[4])
                    for i, v in enumerate(r[6:]):
                        try:
                            p.append(float(v))
                            f.append(lo + i * st)
                        except ValueError:
                            pass
            f, p = np.array(f), np.array(p)
            o = np.argsort(f)
            return f[o], p[o]
        finally:
            os.unlink(path)

    def capture_mpx(self, freq_hz, seconds):
        nbytes = int(MPX_FS * 2 * seconds)
        skip = int(MPX_FS * 2 * 0.3)  # tuner settling
        cmd = [binary("rtl_fm"), "-d", str(self.index), "-M", "fm", "-l", "0", "-A", "std",
               "-p", str(self.ppm), "-s", "171k", "-g", str(self.gain), "-F", "9", "-f", str(int(freq_hz))]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        buf = bytearray()
        deadline = time.time() + seconds + 10
        try:
            while len(buf) < nbytes + skip and time.time() < deadline:
                # A wedged dongle produces no output at all; don't block on read forever.
                if not select.select([p.stdout], [], [], max(deadline - time.time(), 0))[0]:
                    break
                chunk = os.read(p.stdout.fileno(), 65536)
                if not chunk:
                    break
                buf += chunk
        finally:
            p.terminate()
            try:
                p.wait(timeout=3)
            except subprocess.TimeoutExpired:
                p.kill()
        data = bytes(buf[skip:skip + nbytes])
        return np.frombuffer(data[: len(data) // 2 * 2], dtype="<i2")

    def close(self):
        pass


class SoapyBackend:
    """Sweeps by retuning and averaging FFTs in Python, like rtl_power does in C."""
    name = "soapy"

    def __init__(self, args: str = "", gain: float | None = None, ppm: int = 0):
        try:
            import SoapySDR  # noqa: N813
        except ImportError as e:
            raise NoDevice("SoapySDR Python bindings not found (macOS: brew install soapysdr plus the "
                           "module for your radio, e.g. soapyhackrf; Debian/Ubuntu: apt install "
                           "python3-soapysdr soapysdr-module-all; conda: soapysdr)") from e
        self.S = SoapySDR
        found = SoapySDR.Device.enumerate(args)
        if not found:
            raise NoDevice(f"no SoapySDR device matches '{args}'" if args else "no SoapySDR device found")
        self.info = dict(found[0])
        # Open with the same args string: some modules (SoapyRTLSDR) won't match their own
        # enumeration record, or its serial, when it's passed back.
        self.dev = SoapySDR.Device(args)
        rx = SoapySDR.SOAPY_SDR_RX
        self.rate = self._pick_rate()
        self.dev.setSampleRate(rx, 0, self.rate)
        try:
            self.dev.setGainMode(rx, 0, False)  # AGC would change the gain between hops
        except Exception:  # noqa: BLE001
            pass
        rng = self.dev.getGainRange(rx, 0)
        self.gain = gain if gain is not None else round(rng.minimum() + 0.75 * (rng.maximum() - rng.minimum()), 1)
        self.dev.setGain(rx, 0, self.gain)
        if ppm:
            try:
                self.dev.setFrequencyCorrection(rx, 0, ppm)
            except Exception:  # noqa: BLE001
                pass
        self.stream = self.dev.setupStream(rx, SoapySDR.SOAPY_SDR_CF32)
        self.dev.activateStream(self.stream)
        self.buf = np.empty(1 << 16, np.complex64)

    def describe(self) -> str:
        label = self.info.get("label") or self.info.get("driver", "?")
        return f"SoapySDR {label}, {self.rate / 1e6:g} MS/s, gain {self.gain} dB"

    def _pick_rate(self) -> float:
        rx = self.S.SOAPY_SDR_RX
        ranges = self.dev.getSampleRateRange(rx, 0)
        rates = [r.minimum() for r in ranges] + [r.maximum() for r in ranges]

        def ok(x):
            return any(r.minimum() - 1 <= x <= r.maximum() + 1 for r in ranges)
        for want in (2.4e6, 2.5e6, 2.56e6, 3e6, 2.048e6, 2e6):
            if ok(want):
                return want
        above = sorted(r for r in rates if r >= 1.5e6)
        if not above:
            raise NoDevice("device can't sample at 1.5 MS/s or more")
        return min(above[0], 10e6)

    def _read(self, n: int) -> np.ndarray:
        out = np.empty(n, np.complex64)
        got, tries = 0, 0
        while got < n and tries < 200:
            sr = self.dev.readStream(self.stream, [self.buf], min(len(self.buf), n - got), timeoutUs=500_000)
            if sr.ret > 0:
                out[got:got + sr.ret] = self.buf[:sr.ret]
                got += sr.ret
            else:
                tries += 1  # overflow or timeout: keep reading
        return out[:got]

    def _tune(self, f: float) -> None:
        self.dev.setFrequency(self.S.SOAPY_SDR_RX, 0, f)
        self._read(int(self.rate * 0.03))  # discard samples from before/while retuning

    def sweep(self, f_lo, f_hi, bin_hz, seconds, crop=0.3):
        nfft = 1 << int(np.ceil(np.log2(self.rate / bin_hz)))
        df = self.rate / nfft
        keep = int(nfft * (1 - crop)) // 2 * 2
        lo_k = (nfft - keep) // 2
        hop = keep * df
        centres = np.arange(f_lo + hop / 2, f_hi + hop / 2, hop)
        per_hop = max(seconds / len(centres), 16 * nfft / self.rate)
        win = signal.windows.hann(nfft, sym=False).astype(np.float32)
        fs, ps = [], []
        for c in centres:
            self._tune(c)
            x = self._read(int(per_hop * self.rate) // nfft * nfft)
            if len(x) < nfft:
                continue
            x = x[: len(x) // nfft * nfft].reshape(-1, nfft) * win
            P = np.fft.fftshift(np.mean(np.abs(np.fft.fft(x, axis=1)) ** 2, axis=0))
            mid = nfft // 2  # LO leakage (DC spike): replace with its neighbours
            P[mid - 1:mid + 2] = np.median(np.r_[P[mid - 6:mid - 1], P[mid + 2:mid + 7]])
            k = np.arange(lo_k, lo_k + keep)
            fs.append(c + (k - nfft // 2) * df)
            ps.append(10 * np.log10(P[k] + 1e-20))
        if not fs:
            return np.zeros(0), np.zeros(0)
        f, p = np.concatenate(fs), np.concatenate(ps)
        m = (f >= f_lo) & (f <= f_hi)
        o = np.argsort(f[m])
        return f[m][o], p[m][o]

    def capture_mpx(self, freq_hz, seconds):
        offset = 250e3  # tune off the station so the DC spike isn't on the carrier
        self._tune(freq_hz + offset)
        iq = self._read(int(self.rate * seconds))
        t = np.arange(len(iq)) / self.rate
        return fm_demod_mpx(iq * np.exp(2j * np.pi * offset * t).astype(np.complex64), self.rate)

    def close(self):
        try:
            self.dev.deactivateStream(self.stream)
            self.dev.closeStream(self.stream)
        except Exception:  # noqa: BLE001
            pass


def open_backend(spec: str = "auto", gain: float | None = None, ppm: int = 0, log=print):
    """spec: auto | rtl | rtl:INDEX | soapy | soapy:ARGS (e.g. soapy:driver=hackrf)."""
    kind, _, arg = spec.partition(":")
    if kind == "rtl":
        return RtlBackend(gain, ppm, int(arg or 0))
    if kind == "soapy":
        return SoapyBackend(arg, gain, ppm)
    if kind != "auto":
        raise ValueError(f"unknown --sdr '{spec}' (auto, rtl[:INDEX], soapy[:ARGS])")
    errors = []
    if binary("rtl_power") and binary("rtl_fm"):
        name = RtlBackend.probe()
        if name:
            log(f"SDR: {name} via rtl-sdr tools")
            return RtlBackend(gain, ppm)
        errors.append("rtl-sdr tools found no dongle")
    else:
        errors.append("rtl-sdr tools not installed")
    try:
        b = SoapyBackend("", gain, ppm)
        log(f"SDR: {b.describe()}")
        return b
    except NoDevice as e:
        errors.append(str(e))
    raise NoDevice("no SDR found: " + "; ".join(errors))
