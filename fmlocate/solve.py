"""Grid-search position estimate from FM, TV and NOAA Weather Radio channel levels.

Each channel (200 kHz FM, 6 MHz TV, 25 kHz NOAA) is one observation. Its predicted level at a trial position
is the power sum of every licensed transmitter on that frequency (or only the
RDS-identified one), using an empirical field-strength model, plus one unknown
receiver offset `c` (antenna gain, building loss, noise floor) shared by all
channels of one band and integrated out (FM, NOAA, VHF-low TV, VHF-high TV
and UHF TV each get their own, since antenna gain and measurement bandwidth
differ between them). Channels that show nothing are censored
observations: they say "the predicted level here is below the detection
threshold", which is how missing strong stations count against a position.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage
from scipy.special import log_ndtr, logsumexp

from . import fcc, noaa

R_EARTH = 6371.0

# E(dBuV/m) = K0 + 10 log P_kW + KH log HAAT_m - KD log d_km, least-squares fit to the
# FCC 47 CFR 73.211 class reference points (60 dBu contour distance for A, B1, B, C1, C0, C;
# F(50,50) curves). Residuals are under 0.5 dB at those six points.
K0, KH, KD = 89.79, 13.02, 43.91
RX_HEIGHT_M = 9.0
OVER_HORIZON_DB_PER_KM = 0.15  # crude extra loss beyond the radio horizon (a guess, not fitted)


def field_dbu(erp_kw, haat_m, d_km):
    d = np.maximum(d_km, 0.3)
    h = np.clip(haat_m, 30.0, 1000.0)
    lp = 10 * np.log10(np.maximum(erp_kw, 1e-4))
    e = K0 + lp + KH * np.log10(h) - KD * np.log10(d)
    free = 106.92 + lp - 20 * np.log10(d)
    horizon = 4.12 * (np.sqrt(h) + np.sqrt(RX_HEIGHT_M))
    return np.minimum(e, free) - OVER_HORIZON_DB_PER_KM * np.maximum(d - horizon, 0)


def bearing(lat1, lon1, lat2, lon2):
    """Initial great-circle bearing from point 1 to point 2, radians clockwise from north."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(lon2 - lon1)
    return np.arctan2(np.sin(dl) * np.cos(p2), np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl))


def haversine_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = (np.sin((p2 - p1) / 2) ** 2
         + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2)
    return 2 * R_EARTH * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


@dataclass
class Config:
    radius_km: float = 150.0      # search grid half-width around the prior centre
    step_km: float = 2.0
    candidate_margin_km: float = 150.0
    sigma_meas_db: float = 2.0
    sigma_model_db: float = 8.0
    censor_level_db: float = 2.0  # levels below this (above the empty-channel median) = "not received"
    outlier_prob: float = 0.03
    c_range: tuple = (-120.0, -20.0)
    c_step: float = 2.0
    slopes_db: tuple = (-8.0, -4.0, 0.0, 4.0, 8.0)     # per 10 MHz
    direction_amps_db: tuple = (0.0, 8.0, 16.0, 24.0)  # 0 = omnidirectional
    direction_steps: int = 8
    offline: bool = False
    signals: tuple = ("fm", "tv", "wx")


def level_to_snr(level_db):
    """rtl_power level above the floor is ~10 log(1 + S/N); invert it."""
    return 10 * np.log10(np.maximum(10 ** (np.asarray(level_db) / 10) - 1, 1e-3))


def prior_centre(meas: dict, how: str, offline: bool, log=print) -> tuple[float, float, str]:
    if "," in how:
        lat, lon = (float(x) for x in how.split(","))
        return lat, lon, "given on the command line"
    if how in ("auto", "rds"):
        pts = []
        for p in meas.get("rds", []):
            if p.get("callsign"):
                recs = [r for r in fcc.by_call(p["callsign"], offline)
                        if r["service"] == "FM" and abs(r["freq"] - p["freq"]) < 0.05]
                pts += [(r["lat"], r["lon"]) for r in recs]
        if pts:
            lat, lon = np.mean(pts, axis=0)
            return float(lat), float(lon), f"centroid of {len(pts)} RDS-identified transmitters"
        if how == "rds":
            raise RuntimeError("no RDS station identified; use --prior ip or --prior LAT,LON")
    loc = meas.get("ip_location") or ip_location()
    if not loc:
        raise RuntimeError("no RDS station identified and IP geolocation failed; use --prior LAT,LON")
    return loc["lat"], loc["lon"], f"IP geolocation ({loc['location_name']}, coarse guess)"


def ip_location() -> dict | None:
    """Coarse location from the public IP address (ipinfo.io), or None if offline."""
    import json
    import urllib.request
    try:
        req = urllib.request.Request("https://ipinfo.io/json", headers={"User-Agent": fcc.UA})
        with urllib.request.urlopen(req, timeout=8) as r:
            d = json.load(r)
        lat, lon = (float(x) for x in d["loc"].split(","))
        return {"lat": lat, "lon": lon, "location_name": f"{d.get('city')}, {d.get('region')}"}
    except Exception:  # noqa: BLE001
        return None


def _channel_models(meas, cands, extra_by_call, cfg):
    # Most channels are empty, so the median channel level is where "nothing" sits.
    base = float(np.median([c["level_db"] for c in meas["channels"]]))
    levels = {c["freq"]: round(c["level_db"] - base, 2) for c in meas["channels"]}
    rds_call = {p["freq"]: p["callsign"] for p in meas.get("rds", []) if p.get("callsign")}
    by_freq: dict[float, list] = {}
    for r in cands:
        by_freq.setdefault(round(r["freq"], 1), []).append(r)
    chans = []
    for freq, lv in sorted(levels.items()):
        call = rds_call.get(freq)
        neighbour = max(levels.get(round(freq + d, 1), -99) for d in (-0.2, 0.2))
        if not call and neighbour - lv > 10 and lv >= cfg.censor_level_db:
            continue  # probably splatter from the adjacent channel: uninformative
        stations = by_freq.get(freq, [])
        how = "all licences on frequency"
        if call:
            same = [r for r in stations if fcc.base_call(r["call"]) == call]
            same += [r for r in extra_by_call.get(call, [])
                     if round(r["freq"], 1) == freq and r not in same]
            if same:
                stations, how = same, "RDS"
            else:
                relays = [r for r in stations if r["service"] == "FX"]
                if relays:
                    stations, how = relays, "RDS (translator)"
        chans.append({"freq": freq, "label": f"{freq:.1f}", "band": "fm", "level_db": lv,
                      "rds_call": call, "how": how, "stations": stations,
                      "detected": lv >= cfg.censor_level_db,
                      "obs_snr": float(level_to_snr(lv)), "thr_snr": float(level_to_snr(cfg.censor_level_db))})
    return chans


def _tv_models(meas, cands):
    from .measure import tv_channel_levels
    chans = []
    for c in tv_channel_levels(meas["tv"]["bands"]):
        if c["skip"]:
            continue
        det = c["snr_db"] is not None
        chans.append({"freq": c["freq"], "label": f"ch{c['channel']}", "band": c["band"],
                      "level_db": c["level_db"], "rds_call": None,
                      "how": c["how"] or "all licences on channel", "pilot_db": c["pilot_db"],
                      "stations": [r for r in cands if r["channel"] == c["channel"]], "detected": det,
                      "obs_snr": c["snr_db"] if det else c["thr_db"], "thr_snr": c["thr_db"]})
    return chans


def _wx_models(meas, cands, cfg):
    from .measure import wx_channel_levels
    chans = []
    for c in wx_channel_levels(meas["wx"]["spectrum"]):
        lv = c["level_db"]
        chans.append({"freq": c["freq"], "label": f"{c['freq']:.3f}", "band": "wx", "level_db": lv,
                      "rds_call": None, "how": "all NOAA transmitters on frequency",
                      "stations": [r for r in cands if abs(r["freq"] - c["freq"]) < 0.005],
                      "detected": lv >= cfg.censor_level_db,
                      "obs_snr": float(level_to_snr(lv)), "thr_snr": float(level_to_snr(cfg.censor_level_db))})
    return chans


BANDS = ("fm", "wx", "vlo", "vhi", "uhf")


def _combos(cfg) -> list[tuple]:
    combos = [(sl, 0.0, 0.0) for sl in cfg.slopes_db if 0 in cfg.direction_amps_db]
    combos += [(sl, a, phi) for sl in cfg.slopes_db for a in cfg.direction_amps_db if a
               for phi in np.radians(np.arange(0, 360, 360 / cfg.direction_steps))]
    return combos


def _evaluate(chans, lat, lon, cfg, combos):
    """Log posterior (up to a constant) at the given points, integrating out the receiver terms.

    Returns (logpost, best [c, slope, amp, phi] per point, E, B, u).
    """
    E = np.empty((len(chans), lat.size), dtype=np.float32)
    B = np.zeros_like(E)  # bearing (radians) from each point to the channel's dominant transmitter
    for i, ch in enumerate(chans):
        if not ch["stations"]:
            E[i] = -60.0
            continue
        tot = np.zeros(lat.size)
        best = np.full(lat.size, -np.inf)
        for r in ch["stations"]:
            fe = field_dbu(r["erp_kw"], r["haat_m"] or 30, haversine_km(lat, lon, r["lat"], r["lon"]))
            tot += 10 ** (fe / 10)
            better = fe > best
            best[better] = fe[better]
            B[i, better] = bearing(lat[better], lon[better], r["lat"], r["lon"])
        E[i] = 10 * np.log10(tot)

    det = np.array([c["detected"] for c in chans])
    obs = np.array([c["obs_snr"] for c in chans], dtype=np.float32)
    thr = np.array([c["thr_snr"] for c in chans], dtype=np.float32)
    band = np.array([c["band"] for c in chans])
    sigma = float(np.hypot(cfg.sigma_meas_db, cfg.sigma_model_db))
    eps = cfg.outlier_prob

    # Nuisance terms, all integrated out with flat priors: an offset c per band, a slope
    # across the FM band (short whip + building respond differently at 88 and 108 MHz), and a
    # cosine-shaped building attenuation by bearing (amplitude a dB, least loss towards phi)
    # shared by all bands.
    cs = np.arange(cfg.c_range[0], cfg.c_range[1] + 1e-9, cfg.c_step)
    fm = band == "fm"
    u = np.where(fm, (np.array([c["freq"] for c in chans]) - 98.0) / 10.0, 0.0).astype(np.float32)[:, None]
    groups = [b for b in BANDS if (band == b).any()]
    lpost = np.full(lat.size, -np.inf)
    lbest = np.full(lat.size, -np.inf, dtype=np.float32)
    # Per point: slope, amp, phi, then the best c for each band in `groups`.
    best_par = np.zeros((lat.size, 3 + len(groups)), dtype=np.float32)
    log_ok, log_out_d, log_out_c = np.log1p(-eps), np.log(eps / 40.0), np.log(eps)
    for sl, amp, phi in combos:
        Ea = E + sl * u
        if amp:
            Ea = Ea + amp * (np.cos(B - phi) - 1) / 2  # 0 dB towards phi, -amp dB away from it
        total = np.zeros(lat.size)
        tbest = np.zeros(lat.size, dtype=np.float32)
        cbest = []
        for g in groups:
            m = band == g
            Ed, od = Ea[m & det], obs[m & det][:, None]
            Ec, tc = Ea[m & ~det], thr[m & ~det][:, None]
            # Censored channels that can't reach the threshold for any c contribute ~log(1): skip them.
            reach = (Ec.max(1) + cs[-1]) > tc[:, 0] - 4 * sigma
            Ec, tc = Ec[reach], tc[reach]
            ll = np.empty((len(cs), lat.size), dtype=np.float32)
            for k, c in enumerate(cs):
                z = (od - Ed - c) / sigma
                gl = -0.5 * z * z - np.log(sigma * np.sqrt(2 * np.pi))
                a = np.logaddexp(log_ok + gl, log_out_d).sum(0)
                b = np.logaddexp(log_ok + log_ndtr((tc - Ec - c) / sigma), log_out_c).sum(0)
                ll[k] = a + b
            total += logsumexp(ll, axis=0)
            kb = np.argmax(ll, axis=0)
            tbest += ll[kb, np.arange(lat.size)]
            cbest.append(cs[kb])
        lpost = np.logaddexp(lpost, total)
        upd = tbest > lbest
        lbest[upd] = tbest[upd]
        best_par[upd] = np.stack([np.full(upd.sum(), sl), np.full(upd.sum(), amp),
                                  np.full(upd.sum(), phi)] + [c[upd] for c in cbest], 1)
    return lpost, best_par, E, B, u, groups


def solve(meas: dict, prior: str = "auto", cfg: Config | None = None, log=print) -> dict:
    cfg = cfg or Config()
    lat0, lon0, prior_desc = prior_centre(meas, prior, cfg.offline, log)
    log(f"Prior centre: {lat0:.3f}, {lon0:.3f} ({prior_desc}); searching +-{cfg.radius_km:.0f} km")

    cand_radius = int(cfg.radius_km + cfg.candidate_margin_km)
    cands = fcc.band_near(lat0, lon0, cand_radius, cfg.offline)
    extra = {}
    for p in meas.get("rds", []):
        if p.get("callsign"):
            extra[p["callsign"]] = fcc.by_call(p["callsign"], cfg.offline)
    log(f"{len(cands)} FCC FM transmitters within {cand_radius} km")
    tv_cands = wx_cands = []
    if "tv" in cfg.signals and "tv" in meas:
        tv_cands = fcc.tv_near(lat0, lon0, cand_radius, cfg.offline)
        log(f"{len(tv_cands)} FCC TV transmitters within {cand_radius} km")
    if "wx" in cfg.signals and "wx" in meas:
        wx_cands = noaa.near(lat0, lon0, cand_radius, cfg.offline)
        log(f"{len(wx_cands)} NOAA Weather Radio transmitters within {cand_radius} km")

    # Grid in a local equirectangular frame (x east, y north, km).
    n = int(cfg.radius_km / cfg.step_km)
    ax = np.arange(-n, n + 1) * cfg.step_km
    kx = 111.32 * np.cos(np.radians(lat0))
    lats = lat0 + ax / 110.57
    lons = lon0 + ax / kx
    LAT, LON = np.meshgrid(lats, lons, indexing="ij")

    chans = _channel_models(meas, cands, extra, cfg) if "fm" in cfg.signals else []
    if tv_cands:
        chans += _tv_models(meas, tv_cands)
    if wx_cands:
        chans += _wx_models(meas, wx_cands, cfg)
    if not chans:
        raise RuntimeError("no channels to use (check --signals and the measurement)")
    det = np.array([c["detected"] for c in chans])
    combos = _combos(cfg)
    by_band = ", ".join(f"{b} {sum(c['detected'] for c in chans if c['band'] == b)}/"
                        f"{sum(c['band'] == b for c in chans)}" for b in BANDS if any(c["band"] == b for c in chans))
    log(f"{len(chans)} channels used ({det.sum()} received, {(~det).sum()} censored; {by_band}), "
        f"sigma {np.hypot(cfg.sigma_meas_db, cfg.sigma_model_db):.1f} dB, "
        f"{len(combos)} receiver-response hypotheses")

    # Coarse pass at 3x the step over the whole area, then the fine grid only where the
    # coarse posterior is within e^-12 of its peak (plus one coarse cell of margin).
    m = max(1, int(round(n / 3)))
    cax = np.linspace(-n, n, 2 * m + 1) * cfg.step_km
    CLAT, CLON = np.meshgrid(lat0 + cax / 110.57, lon0 + cax / kx, indexing="ij")
    lc, *_ = _evaluate(chans, CLAT.ravel(), CLON.ravel(), cfg, combos)
    keep_c = (lc > lc.max() - 12).reshape(CLAT.shape)
    keep_c = ndimage.binary_dilation(keep_c, iterations=1)
    ci = np.clip(np.round((ax + n * cfg.step_km) / (cax[1] - cax[0])).astype(int), 0, 2 * m)
    keep = keep_c[np.ix_(ci, ci)].ravel()
    idx = np.flatnonzero(keep)
    log(f"Refining {idx.size} of {LAT.size} cells at {cfg.step_km:g} km")
    lsub, par_sub, E, B, u, groups = _evaluate(chans, LAT.ravel()[idx], LON.ravel()[idx], cfg, combos)
    logpost = np.full(LAT.size, -np.inf)
    logpost[idx] = lsub
    sigma = float(np.hypot(cfg.sigma_meas_db, cfg.sigma_model_db))

    w = np.exp(logpost - logpost.max())
    w /= w.sum()
    X = np.broadcast_to(ax[None, :], LAT.shape).ravel()
    Y = np.broadcast_to(ax[:, None], LAT.shape).ravel()
    mx, my = float(w @ X), float(w @ Y)
    cov = np.cov(np.vstack([X, Y]), aweights=w)
    evals, evecs = np.linalg.eigh(cov)
    i_map = int(np.argmax(w))
    order = np.argsort(-w)
    cum = np.cumsum(w[order])
    cell = cfg.step_km ** 2

    def hpd_radius(q):
        return float(np.sqrt((np.searchsorted(cum, q) + 1) * cell / np.pi))

    hpd90 = np.zeros(w.size, bool)
    hpd90[order[: np.searchsorted(cum, 0.9) + 1]] = True
    near_edge = float(w[(np.abs(X) > cfg.radius_km - 3 * cfg.step_km)
                        | (np.abs(Y) > cfg.radius_km - 3 * cfg.step_km)].sum())

    est = {"lat": float(LAT.ravel()[i_map]), "lon": float(LON.ravel()[i_map])}
    mean = {"lat": lat0 + my / 110.57, "lon": lon0 + mx / kx}
    j_map = int(np.searchsorted(idx, i_map))  # position within the refined subset
    slope, amp, phi = (float(x) for x in par_sub[j_map][:3])
    offsets = {g: float(x) for g, x in zip(groups, par_sub[j_map][3:])}

    table = []
    for i, ch in enumerate(chans):
        pred_e = float(E[i, j_map]) + slope * float(u[i, 0])
        pred_e += amp * (np.cos(float(B[i, j_map]) - phi) - 1) / 2
        pred_snr = pred_e + offsets[ch["band"]]
        if not ch["detected"] and pred_snr < ch["thr_snr"] - 6:
            continue
        st = None
        if ch["stations"]:
            d = [haversine_km(est["lat"], est["lon"], r["lat"], r["lon"]) for r in ch["stations"]]
            fe = [field_dbu(r["erp_kw"], r["haat_m"] or 30, dd) for r, dd in zip(ch["stations"], d)]
            j = int(np.argmax(fe))
            r = ch["stations"][j]
            st = {"call": r["call"], "service": r["service"], "city": f"{r['city']}, {r['state']}",
                  "lat": r["lat"], "lon": r["lon"], "erp_kw": r["erp_kw"], "haat_m": r["haat_m"],
                  "distance_km": round(d[j], 1)}
        table.append({"freq": ch["freq"], "label": ch["label"], "band": ch["band"],
                      "level_db": ch["level_db"], "rds_call": ch["rds_call"],
                      "identified_by": ch["how"], "n_licences": len(ch["stations"]),
                      "detected": ch["detected"], "obs_snr_db": round(ch["obs_snr"], 1) if ch["detected"] else None,
                      "pred_snr_db": round(pred_snr, 1), "station": st})

    return {
        "prior": {"lat": lat0, "lon": lon0, "source": prior_desc, "radius_km": cfg.radius_km},
        "estimate": est, "posterior_mean": mean,
        "ellipse_1sigma_km": {"major": float(np.sqrt(evals[1])), "minor": float(np.sqrt(evals[0])),
                              "major_bearing_deg": float(np.degrees(np.arctan2(evecs[0, 1], evecs[1, 1])) % 180)},
        "radius50_km": hpd_radius(0.5), "radius90_km": hpd_radius(0.9),
        "mass_near_grid_edge": near_edge,
        "receiver_offset_db": offsets.get("fm"),
        "band_offsets_db": offsets,
        "receiver_response": {"slope_db_per_10mhz": slope, "direction_amp_db": amp,
                              "least_loss_bearing_deg": float(np.degrees(phi)) % 360}, "sigma_db": sigma, "config": cfg.__dict__,
        "channels": table,
        "grid": {"lats": lats.tolist(), "lons": lons.tolist(), "post": w.reshape(LAT.shape),
                 "hpd90": hpd90.reshape(LAT.shape)},
        "transmitters": [{k: r[k] for k in ("call", "freq", "service", "lat", "lon", "erp_kw", "haat_m")}
                         for ch in chans if ch["detected"] for r in ch["stations"]],
    }


def in_hpd90(res: dict, lat: float, lon: float) -> bool:
    g = res["grid"]
    i = int(np.argmin(np.abs(np.array(g["lats"]) - lat)))
    j = int(np.argmin(np.abs(np.array(g["lons"]) - lon)))
    inside = abs(g["lats"][i] - lat) < 0.05 and abs(g["lons"][j] - lon) < 0.05
    return bool(inside and g["hpd90"][i, j])
