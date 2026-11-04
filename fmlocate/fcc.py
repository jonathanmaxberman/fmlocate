"""FCC FM and TV licence data (transmitter site, ERP, HAAT), with an on-disk cache.

FM uses the fmq endpoint, parsing ERP and HAAT, and fetches the whole band around a point in one request. TV uses tvq, which
returns the same column layout with the RF channel in place of the frequency.
"""
from __future__ import annotations

import hashlib
import re
import subprocess
import time
import urllib.parse

from . import paths

URL = "https://transition.fcc.gov/fcc-bin/fmq?"
TV_URL = "https://transition.fcc.gov/fcc-bin/tvq?"
UA = "fmlocate/0.2 (https://github.com/jonathanmaxberman/fmlocate)"
# Full-service, translators, low-power and boosters. FS/FA are auxiliary/application records.
SERVICES = {"FM", "FX", "FL", "FB"}
# Full-service digital, distributed-transmission sites, Class A, low-power and translators.
TV_SERVICES = {"DTV", "DTS", "DCA", "LPD", "DRT"}


def _get(params: dict, offline: bool = False, url: str = URL) -> str:
    q = urllib.parse.urlencode(params)
    key = q if url == URL else url + q
    path = paths.cache("fcc") / (hashlib.sha1(key.encode()).hexdigest()[:16] + ".txt")
    if path.exists():
        return path.read_text()
    if offline:
        raise FileNotFoundError(f"FCC query not cached and --offline given: {q}")
    # The FCC front end rejects some urllib requests that curl gets through with, so use curl.
    for attempt in range(4):
        p = subprocess.run(["curl", "-sSf", "--max-time", "90", "-A", UA, url + q],
                           capture_output=True, text=True)
        # The FCC server answers 503 now and then, especially to several queries at once.
        if p.returncode == 0 or "error: 5" not in p.stderr:
            break
        time.sleep(10 * (attempt + 1))
    if p.returncode:
        raise RuntimeError(f"FCC query failed: {p.stderr.strip()}")
    text = p.stdout
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return text


def _num(s: str) -> float | None:
    m = re.match(r"\s*(-?[\d.]+)", s)
    return float(m.group(1)) if m else None


def parse(text: str) -> list[dict]:
    out = []
    for line in text.splitlines():
        f = [x.strip() for x in line.split("|")]
        if len(f) < 28 or not f[1]:
            continue
        try:
            lat = float(f[20]) + float(f[21]) / 60 + float(f[22]) / 3600
            lon = float(f[24]) + float(f[25]) / 60 + float(f[26]) / 3600
        except ValueError:
            continue
        erp = [x for x in (_num(f[14]), _num(f[15])) if x is not None]
        haat = [x for x in (_num(f[16]), _num(f[17])) if x is not None]
        out.append({
            "call": f[1], "freq": _num(f[2]), "service": f[3], "directional": f[5] == "DA",
            "status": f[9], "city": f[10].title(), "state": f[11], "facility": f[18],
            "erp_kw": max(erp) if erp else None, "haat_m": max(haat) if haat else None,
            "lat": -lat if f[19] == "S" else lat, "lon": -lon if f[23] == "W" else lon,
        })
    return out


def _dms(v: float) -> tuple[int, int, int]:
    t = int(round(abs(v) * 3600))  # whole seconds, so 39.8 gives 39 48 0, not 39 47 60
    return t // 3600, t // 60 % 60, t % 60


def _best_per_facility(recs: list[dict]) -> list[dict]:
    rank = {"LIC": 0, "STA": 1, "MOD": 2}
    best: dict[tuple, dict] = {}
    for r in recs:
        if r["service"] not in SERVICES or r["status"] not in rank or not r["erp_kw"]:
            continue
        key = (r["facility"] or r["call"], r["freq"])
        if key not in best or rank[r["status"]] < rank[best[key]["status"]]:
            best[key] = r
    return list(best.values())


def band_near(lat: float, lon: float, radius_km: int, offline: bool = False) -> list[dict]:
    """Every licensed FM transmitter within radius_km of a point."""
    d1, m1, s1 = _dms(lat)
    d2, m2, s2 = _dms(lon)
    text = _get({"freq": "87.9", "fre2": "107.9", "list": "4", "size": "9", "dist": str(radius_km),
                 "NS": "N" if lat >= 0 else "S", "dlat2": d1, "mlat2": m1, "slat2": s1,
                 "EW": "W" if lon < 0 else "E", "dlon2": d2, "mlon2": m2, "slon2": s2}, offline)
    return _best_per_facility(parse(text))


def tv_channel_low_mhz(ch: int) -> float:
    """Lower edge of a US 6 MHz TV channel (2-36 after the 600 MHz repack)."""
    if ch <= 4:
        return 54.0 + 6 * (ch - 2)
    if ch <= 6:
        return 76.0 + 6 * (ch - 5)
    if ch <= 13:
        return 174.0 + 6 * (ch - 7)
    return 470.0 + 6 * (ch - 14)


TV_CHANNELS = list(range(2, 37))


def tv_near(lat: float, lon: float, radius_km: int, offline: bool = False) -> list[dict]:
    """Every licensed digital TV transmitter within radius_km of a point.

    `freq` is the channel centre in MHz. Licences sharing one transmitter (channel
    sharing, or a DTS site listed again as the main DTV record) are one signal, so they
    are merged per channel and site.
    """
    d1, m1, s1 = _dms(lat)
    d2, m2, s2 = _dms(lon)
    text = _get({"chan": "2", "cha2": "36", "list": "4", "size": "9", "dist": str(radius_km),
                 "NS": "N" if lat >= 0 else "S", "dlat2": d1, "mlat2": m1, "slat2": s1,
                 "EW": "W" if lon < 0 else "E", "dlon2": d2, "mlon2": m2, "slon2": s2},
                offline, TV_URL)
    rank = {"LIC": 0, "STA": 1}
    best: dict[tuple, dict] = {}
    for line in text.splitlines():
        f = [x.strip() for x in line.split("|")]
        if len(f) < 37 or not f[1] or f[3] not in TV_SERVICES or f[9] not in rank:
            continue
        try:
            ch = int(f[4])
            lat_ = float(f[20]) + float(f[21]) / 60 + float(f[22]) / 3600
            lon_ = float(f[24]) + float(f[25]) / 60 + float(f[26]) / 3600
        except ValueError:
            continue
        erp = _num(f[14])
        # Some Class A / low-power records leave column 16 at 0; column 36 is the HAAT there.
        haat = _num(f[16]) or _num(f[36])
        if not erp or not 2 <= ch <= 36:
            continue
        r = {"call": f[1], "channel": ch, "freq": tv_channel_low_mhz(ch) + 3.0, "service": f[3],
             "directional": f[5] == "DA", "status": f[9], "city": f[10].title(), "state": f[11],
             "facility": f[18], "erp_kw": erp, "haat_m": haat,
             "lat": -lat_ if f[19] == "S" else lat_, "lon": -lon_ if f[23] == "W" else lon_}
        key = (ch, round(r["lat"], 3), round(r["lon"], 3))
        old = best.get(key)
        if old is None or (rank[r["status"]], -erp) < (rank[old["status"]], -old["erp_kw"]):
            best[key] = r
    return list(best.values())


def by_call(call: str, offline: bool = False) -> list[dict]:
    return _best_per_facility(parse(_get({"call": call, "list": "4", "size": "9"}, offline)))


def base_call(call: str) -> str:
    return re.split(r"[-\s]", call.strip().upper())[0]
