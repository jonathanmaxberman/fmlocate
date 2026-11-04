"""Heatmap PNG of the position estimate, on an OpenStreetMap background when online."""
from __future__ import annotations

import io
import math
import urllib.request

import numpy as np

from . import paths
from .fcc import UA

TILE_URL = "https://tile.openstreetmap.org/{z}/{x}/{y}.png"
R = 6378137.0


def merc(lat, lon):
    """Web Mercator metres."""
    lat = np.clip(np.asarray(lat, float), -85, 85)
    return R * np.radians(lon), R * np.log(np.tan(np.pi / 4 + np.radians(lat) / 2))


def _tile(z: int, x: int, y: int):
    import matplotlib.image as mpimg
    path = paths.cache("tiles", str(z), str(x), f"{y}.png")
    if not path.exists():
        req = urllib.request.Request(TILE_URL.format(z=z, x=x, y=y), headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=10) as r:
            data = r.read()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    return mpimg.imread(io.BytesIO(path.read_bytes()), format="png")


def _basemap(ax, x0, x1, y0, y1, width_px: int) -> bool:
    """Draw greyed OSM tiles covering the extent. False if they couldn't be fetched."""
    span = x1 - x0
    z = int(np.clip(math.floor(math.log2(2 * math.pi * R * width_px / (256 * span))), 1, 14))
    n = 2 ** z

    def tx(x):
        return (x + math.pi * R) / (2 * math.pi * R) * n

    def ty(y):
        return (math.pi * R - y) / (2 * math.pi * R) * n
    xs = range(int(tx(x0)), int(tx(x1)) + 1)
    ys = range(int(ty(y1)), int(ty(y0)) + 1)
    if len(xs) * len(ys) > 64:
        return False
    try:
        rows = []
        for y in ys:
            row = []
            for x in xs:
                t = _tile(z, x % n, y)
                row.append(t[..., :3] if t.ndim == 3 else np.dstack([t] * 3))
            rows.append(np.concatenate(row, axis=1))
        img = np.concatenate(rows, axis=0)
    except Exception:  # noqa: BLE001  (offline, blocked, bad tile)
        return False
    grey = img[..., :3] @ np.array([0.3, 0.59, 0.11])
    grey = 0.55 + 0.45 * grey  # lighten so the heat stays readable
    size = 2 * math.pi * R / n
    ext = (xs[0] * size - math.pi * R, (xs[-1] + 1) * size - math.pi * R,
           math.pi * R - (ys[-1] + 1) * size, math.pi * R - ys[0] * size)
    ax.imshow(grey, extent=ext, cmap="gray", vmin=0, vmax=1, zorder=0, interpolation="bilinear")
    return True


def _hpd_levels(w: np.ndarray, qs=(0.5, 0.9)) -> list[float]:
    v = np.sort(w.ravel())[::-1]
    c = np.cumsum(v)
    return [float(v[min(np.searchsorted(c, q), len(v) - 1)]) for q in qs]


def write_png(res: dict, meas: dict, path: str, truth: dict | None = None, basemap: bool = True,
              title: str | None = None) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    g = res["grid"]
    lats, lons, w = np.array(g["lats"]), np.array(g["lons"]), np.asarray(g["post"])
    LAT, LON = np.meshgrid(lats, lons, indexing="ij")

    # Extent: the cells holding 99% of the probability, padded, at least 60 km across.
    lev99 = _hpd_levels(w, (0.99,))[0]
    m = w >= lev99
    pts_lat, pts_lon = list(LAT[m]), list(LON[m])
    if truth:
        pts_lat.append(truth["lat"])
        pts_lon.append(truth["lon"])
    X, Y = merc(pts_lat, pts_lon)
    cx, cy = (X.min() + X.max()) / 2, (Y.min() + Y.max()) / 2
    half = max(X.max() - X.min(), Y.max() - Y.min()) * 0.75
    half = max(half, 30e3 / math.cos(math.radians(res["estimate"]["lat"])))
    x0, x1, y0, y1 = cx - half * 1.25, cx + half * 1.25, cy - half, cy + half

    fig, ax = plt.subplots(figsize=(10, 8), dpi=120)
    have_map = basemap and _basemap(ax, x0, x1, y0, y1, 1200)
    if not have_map:
        ax.set_facecolor("#f4f4f2")

    # Grid cell edges in Mercator (the grid is regular in degrees).
    dlat, dlon = lats[1] - lats[0], lons[1] - lons[0]
    elat = np.r_[lats - dlat / 2, lats[-1] + dlat / 2]
    elon = np.r_[lons - dlon / 2, lons[-1] + dlon / 2]
    EX, EY = merc(*np.meshgrid(elat, elon, indexing="ij"))
    cmap = LinearSegmentedColormap.from_list(
        "heat", [(1, 0.95, 0.4, 0), (1, 0.75, 0.1, 0.45), (0.95, 0.35, 0.05, 0.7), (0.65, 0, 0.1, 0.85)])
    v = np.sqrt(w / w.max())  # square root keeps the tails visible
    ax.pcolormesh(EX, EY, np.ma.masked_less(v, 0.02), cmap=cmap, vmin=0, vmax=1, shading="flat", zorder=2)
    CX, CY = merc(LAT, LON)
    l50, l90 = _hpd_levels(w)
    cs = ax.contour(CX, CY, w, levels=sorted([l90, l50]), colors=["#7a0019", "#7a0019"],
                    linewidths=[1.2, 2.0], linestyles=["--", "-"], zorder=3)
    ax.clabel(cs, fmt={min(l90, l50): "90%", max(l90, l50): "50%"}, fontsize=8)

    # Transmitters on received channels, grouped by site.
    from .solve import haversine_km
    sites: dict[tuple, list] = {}
    for c in res["channels"]:
        s = c["station"]
        if c["detected"] and s:
            key = next((k for k in sites if haversine_km(k[0], k[1], s["lat"], s["lon"]) < 3),
                       (s["lat"], s["lon"]))
            sites.setdefault(key, []).append((c, s))
    for (la, lo), items in sites.items():
        sx, sy = merc(la, lo)
        if not (x0 < sx < x1 and y0 < sy < y1):
            continue
        ax.plot(sx, sy, "^", ms=7, mfc="#1565c0", mec="white", mew=0.8, zorder=5)
        calls = sorted({s["call"] for _, s in items})
        label = calls[0] + (f" +{len(calls) - 1}" if len(calls) > 1 else "")
        ax.annotate(label, (sx, sy), xytext=(5, 4), textcoords="offset points", fontsize=7,
                    color="#0d3c75", zorder=6)

    e, pm = res["estimate"], res["posterior_mean"]
    ex, ey = merc(e["lat"], e["lon"])
    ax.plot(ex, ey, "X", ms=12, mfc="#c62828", mec="white", mew=1.2, zorder=7, label="Most likely")
    mx, my = merc(pm["lat"], pm["lon"])
    ax.plot(mx, my, "P", ms=9, mfc="#6a1b9a", mec="white", mew=1, zorder=7, label="Posterior mean")
    if truth:
        tx_, ty_ = merc(truth["lat"], truth["lon"])
        ax.plot(tx_, ty_, "o", ms=13, mfc="none", mec="#2e7d32", mew=2.5, zorder=7, label="Known location")
    ax.plot([], [], "^", ms=7, mfc="#1565c0", mec="white", label="Received transmitter site")
    ax.plot([], [], color="#7a0019", lw=2, label="50% / 90% region")

    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    # Degree ticks, so the image is usable without a basemap too.
    lon0, lon1 = np.degrees(x0 / R), np.degrees(x1 / R)
    lat0, lat1 = (np.degrees(2 * np.arctan(np.exp(y / R)) - np.pi / 2) for y in (y0, y1))
    st = min((0.1, 0.2, 0.25, 0.5, 1, 2), key=lambda d: abs((lon1 - lon0) / d - 5))
    tlon = np.arange(np.ceil(lon0 / st) * st, lon1, st)
    tlat = np.arange(np.ceil(lat0 / st) * st, lat1, st)
    ax.set_xticks(merc(0, tlon)[0], [f"{abs(v):.2f}°{'W' if v < 0 else 'E'}" for v in tlon], fontsize=7)
    ax.set_yticks(merc(tlat, 0)[1], [f"{abs(v):.2f}°{'S' if v < 0 else 'N'}" for v in tlat], fontsize=7)
    ax.grid(color="#888", lw=0.4, alpha=0.4, zorder=1)
    ax.legend(loc="lower left", fontsize=8, framealpha=0.9)

    # Scale bar (Mercator metres are stretched by 1/cos(lat)).
    k = 1 / math.cos(math.radians(e["lat"]))
    km = 10 if half / k < 40e3 else 25 if half / k < 100e3 else 50
    bx, by = x1 - (x1 - x0) * 0.05 - km * 1e3 * k, y0 + (y1 - y0) * 0.04
    ax.plot([bx, bx + km * 1e3 * k], [by, by], color="k", lw=3, zorder=8)
    ax.text(bx + km * 500 * k, by + (y1 - y0) * 0.012, f"{km} km", ha="center", fontsize=8, zorder=8)

    bands = sorted({c["band"] for c in res["channels"]}, key=["fm", "wx", "vlo", "vhi", "uhf"].index)
    names = {"fm": "FM", "wx": "NOAA", "vlo": "TV", "vhi": "TV", "uhf": "TV"}
    used = "+".join(dict.fromkeys(names[b] for b in bands))
    lines = [f"Most likely  {e['lat']:.4f}, {e['lon']:.4f}",
             f"Mean         {pm['lat']:.4f}, {pm['lon']:.4f}",
             f"50% within ~{res['radius50_km']:.0f} km, 90% within ~{res['radius90_km']:.0f} km",
             f"Signals: {used}"]
    if truth:
        lines.append(f"Known location {truth['err']:.1f} km from most likely")
    ax.text(0.01, 0.99, "\n".join(lines), transform=ax.transAxes, va="top", ha="left", family="monospace",
            fontsize=8.5, bbox=dict(boxstyle="round", fc="white", ec="#999", alpha=0.92), zorder=9)
    when = meas.get("time_utc", "")
    ax.set_title(title or f"fmlocate estimate  {when}", fontsize=11)
    credit = "Map © OpenStreetMap contributors" if have_map else "(no basemap)"
    ax.text(0.99, 0.005, credit, transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5,
            color="#333", zorder=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
