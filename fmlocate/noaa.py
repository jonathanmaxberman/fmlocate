"""NOAA Weather Radio transmitters (site, frequency, power), with an on-disk cache.

The NWS station pages load one nationwide list, ccl-data.js. It gives each
transmitter's coordinates and power but not its antenna height, so a typical
height is assumed.
"""
from __future__ import annotations

import json
import subprocess
from . import paths
from .fcc import UA

URL = "https://www.weather.gov/source/nwr/JS/ccl-data.js"
CHANNELS_MHZ = [162.400, 162.425, 162.450, 162.475, 162.500, 162.525, 162.550]
ASSUMED_HAAT_M = 150.0  # not published; NWR antennas are mostly on 100-300 m broadcast towers


def _load(offline: bool) -> list[dict]:
    CACHE = paths.cache("noaa", "ccl-data.json")
    if CACHE.exists():
        return json.loads(CACHE.read_text())
    if offline:
        raise FileNotFoundError("NOAA transmitter list not cached and --offline given")
    p = subprocess.run(["curl", "-sSf", "--max-time", "60", "-A", UA, URL], capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError(f"NOAA station list download failed: {p.stderr.strip()}")
    t = p.stdout
    data = json.loads(t[t.index("["):t.rindex("]") + 1])
    out = []
    for x in data:
        try:
            out.append({"call": x["callsign"], "freq": round(float(x["freq"]), 3), "service": "NWR",
                        "status": x["status"], "city": x.get("sitename") or x.get("siteloc", ""),
                        "state": x.get("sitestate", ""), "erp_kw": float(x["power"]) / 1000,
                        "haat_m": ASSUMED_HAAT_M, "lat": float(x["lat"]), "lon": float(x["lon"])})
        except (KeyError, ValueError):
            continue
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(out))
    return out


def near(lat: float, lon: float, radius_km: float, offline: bool = False) -> list[dict]:
    """Transmitters in service within radius_km. Power is transmitter output, used as ERP."""
    from .solve import haversine_km
    return [r for r in _load(offline) if r["status"] != "OUT OF SERVICE" and r["erp_kw"] > 0
            and haversine_km(lat, lon, r["lat"], r["lon"]) <= radius_km]
