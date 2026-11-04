"""fmlocate: estimate where this computer is from the FM, TV and NOAA Weather Radio stations an
attached SDR receives. No GPS, Wi-Fi or internet location service needed (the FCC and NOAA
transmitter lists are downloaded once and cached)."""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

from . import __version__, paths, solve
from .solve import Config, haversine_km


def _latlon(s: str) -> dict:
    lat, lon = (float(x) for x in s.split(","))
    return {"lat": lat, "lon": lon}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="fmlocate", description=__doc__)
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--replay", metavar="JSON", help="use a saved measurement instead of the SDR")
    src.add_argument("--live", action="store_true", help="if the tuner is busy (RadioShelf), fail instead "
                     "of replaying the newest saved measurement")
    src.add_argument("--list-sdrs", action="store_true", help="list the SDRs fmlocate can see, then exit")
    ap.add_argument("--sdr", default="auto", help="receiver: auto, rtl[:INDEX] (rtl-sdr tools), or "
                    "soapy[:ARGS] (any SoapySDR device, e.g. soapy:driver=hackrf)  [auto]")
    ap.add_argument("--ppm", type=int, default=0, help="tuner frequency correction, ppm [0]")
    ap.add_argument("--save", metavar="JSON", help="where to save a live measurement "
                    f"(default {paths.measurements()}/<time>.json)")
    ap.add_argument("--prior", default="auto", help="search centre: auto (RDS transmitters, else IP), "
                    "rds, ip, or LAT,LON  [auto]")
    ap.add_argument("--radius", type=float, default=150, help="search half-width, km [150]")
    ap.add_argument("--step", type=float, default=2, help="grid step, km [2]")
    ap.add_argument("--sigma-model", type=float, default=8, help="propagation + building error, dB [8]")
    ap.add_argument("--simple", action="store_true", help="only fit a common receiver offset "
                    "(no frequency slope or directional building loss)")
    ap.add_argument("--gain", type=float, help="tuner gain, dB [rtl: 40.2; soapy: 75%% of the device's range]")
    ap.add_argument("--signals", default="fm,tv,wx", help="which signals to use: any of fm, tv, wx "
                    "(NOAA Weather Radio). FM is always swept, for RDS and the search centre [fm,tv,wx]")
    ap.add_argument("--offline", action="store_true", help="use only cached FCC/NOAA data (and no basemap)")
    ap.add_argument("--cache-dir", metavar="DIR", help=f"licence and map-tile cache [{paths.CACHE_DIR}]")
    ap.add_argument("--truth", metavar="LAT,LON", help="known location, to report the error")
    ap.add_argument("--png", metavar="FILE", help="heatmap image path [fmlocate-<time>.png]")
    ap.add_argument("--no-png", action="store_true", help="don't write the heatmap image")
    ap.add_argument("--no-basemap", action="store_true", help="heatmap without OpenStreetMap tiles")
    ap.add_argument("--map", metavar="HTML", help="also write an interactive Leaflet map")
    ap.add_argument("--json", metavar="OUT", help="write the result as JSON")
    ap.add_argument("-q", "--quiet", action="store_true", help="print only the report")
    ap.add_argument("--version", action="version", version=f"fmlocate {__version__}")
    a = ap.parse_args(argv)
    if a.cache_dir:
        paths.set_cache_dir(a.cache_dir)
    signals = tuple(x.strip() for x in a.signals.split(",") if x.strip())
    if set(signals) - {"fm", "tv", "wx"}:
        ap.error("--signals takes fm, tv and/or wx")
    log = (lambda *x, **k: None) if a.quiet else (lambda *x, **k: print(*x, **k, file=sys.stderr, flush=True))

    if a.list_sdrs:
        return list_sdrs()
    if a.replay:
        meas = json.loads(Path(a.replay).read_text())
        log(f"Replaying {a.replay} ({meas.get('time_utc')})")
    else:
        from .measure import TunerBusy, measure
        from .sdr import NoDevice
        try:
            meas = measure(a.sdr, gain=a.gain, ppm=a.ppm, extra=tuple(x for x in signals if x != "fm"), log=log)
        except NoDevice as e:
            log(f"{e}\nUse --list-sdrs to check, or --replay FILE to use a saved measurement.")
            return 2
        except TunerBusy as e:
            saved = sorted(paths.measurements().glob("*.json"))
            if a.live or not saved:
                log(f"Tuner not available: {e}")
                return 2
            log(f"Tuner not available ({e}); replaying newest saved sweep {saved[-1].name}")
            meas = json.loads(saved[-1].read_text())
        else:
            out = Path(a.save) if a.save else paths.measurements() / (
                datetime.now().strftime("%Y%m%d-%H%M%S") + ".json")
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(json.dumps(meas, indent=1))
            log(f"Saved measurements to {out}")

    if "spectrum" in meas:
        from .measure import levels_from_spectrum
        meas["channels"] = levels_from_spectrum(meas)  # current level definition, not the saved one
    cfg = Config(radius_km=a.radius, step_km=a.step, sigma_model_db=a.sigma_model, offline=a.offline,
                 signals=signals)
    missing = [x for x in signals if x != "fm" and x not in meas]
    if missing:
        log(f"(this measurement has no {', '.join(missing)} data)")
    if a.simple:
        cfg.slopes_db, cfg.direction_amps_db = (0.0,), (0.0,)
    res = solve.solve(meas, a.prior, cfg, log)
    print(report(res, a.truth))
    truth = None
    if a.truth:
        t = _latlon(a.truth)
        truth = {**t, "err": float(haversine_km(res["estimate"]["lat"], res["estimate"]["lon"], t["lat"], t["lon"]))}
    if not a.no_png:
        from . import plot
        png = a.png or f"fmlocate-{datetime.now().strftime('%Y%m%d-%H%M%S')}.png"
        plot.write_png(res, meas, png, truth, basemap=not (a.no_basemap or a.offline))
        log(f"Heatmap written to {png}")
    if a.map:
        from . import mapview
        mapview.write(res, meas, a.map, truth)
        log(f"Map written to {a.map}")
    if a.json:
        out = {k: v for k, v in res.items() if k != "grid"}
        Path(a.json).write_text(json.dumps(out, indent=1))
    return 0


BAND_NAMES = {"fm": "FM", "wx": "NOAA", "vlo": "TV VHF-low", "vhi": "TV VHF-high", "uhf": "TV UHF"}


def list_sdrs() -> int:
    from .sdr import RtlBackend, binary
    print("rtl-sdr tools:", "installed" if binary("rtl_power") else "not installed")
    if binary("rtl_test"):
        print("  device 0:", RtlBackend.probe() or "none found")
    try:
        import SoapySDR  # noqa: N813
    except ImportError:
        print("SoapySDR Python bindings: not installed")
        return 0
    devs = SoapySDR.Device.enumerate()
    print(f"SoapySDR: {len(devs)} device(s)")
    for d in devs:
        d = dict(d)
        args = ",".join(f"{k}={v}" for k, v in d.items() if k in ("driver", "serial"))
        print(f"  {d.get('label', d.get('driver'))}   --sdr soapy:{args}")
    return 0


def report(res: dict, truth: str | None) -> str:
    e, m, el = res["estimate"], res["posterior_mean"], res["ellipse_1sigma_km"]
    lines = [
        f"Estimate:        {e['lat']:.4f}, {e['lon']:.4f}   (most likely grid cell)",
        f"Posterior mean:  {m['lat']:.4f}, {m['lon']:.4f}",
        f"Uncertainty:     50% within ~{res['radius50_km']:.0f} km, 90% within ~{res['radius90_km']:.0f} km "
        f"(area-equivalent radii)",
        f"1-sigma ellipse: {el['major']:.0f} x {el['minor']:.0f} km, long axis bearing {el['major_bearing_deg']:.0f} deg",
        f"Search prior:    {res['prior']['source']}, +-{res['prior']['radius_km']:.0f} km",
        "Receiver offset: " + ", ".join(f"{BAND_NAMES[b]} {v:.0f}" for b, v in res["band_offsets_db"].items())
        + " dB (S/N = field dBuV/m + offset)",
        (("Receiver fit:    {slope_db_per_10mhz:+.0f} dB per 10 MHz across FM; " if "fm" in res["band_offsets_db"]
          else "Receiver fit:    ") +
         "building loss up to {direction_amp_db:.0f} dB, least towards {least_loss_bearing_deg:.0f} deg")
        .format(**res["receiver_response"]),
    ]
    if res["mass_near_grid_edge"] > 0.05:
        lines.append(f"WARNING: {100 * res['mass_near_grid_edge']:.0f}% of the probability is at the edge "
                     "of the search area; widen --radius or change --prior.")
    if truth:
        t = _latlon(truth)
        d1 = haversine_km(e["lat"], e["lon"], t["lat"], t["lon"])
        d2 = haversine_km(m["lat"], m["lon"], t["lat"], t["lon"])
        lines.append(f"Known location:  estimate is {d1:.1f} km away (posterior mean {d2:.1f} km); "
                     f"{'inside' if solve.in_hpd90(res, t['lat'], t['lon']) else 'OUTSIDE'} the 90% region")
    lines += ["", "channel level   S/N obs  pred   station (dominant licence at estimate)          dist"]
    band = None
    for c in res["channels"]:
        if c["band"] != band and len(res["band_offsets_db"]) > 1:
            band = c["band"]
            lines.append(f"  -- {BAND_NAMES[band]}")
        s = c["station"]
        name = (f"{s['call']:<10} {s['erp_kw']:>6.2f} kW {s['haat_m'] or 0:>4.0f} m  {s['city'][:18]:<18}"
                if s else "(no licence)")
        tag = "RDS" if c["rds_call"] else f"{c['n_licences']:>2}x "
        obs = f"{c['obs_snr_db']:>5.1f}" if c["detected"] else "    -"
        dist = f"{s['distance_km']:>5.0f} km" if s else ""
        lines.append(f"{c['label']:>7} {c['level_db']:>5.1f}  {obs}  {c['pred_snr_db']:>5.1f}   {tag} {name} {dist}")
    lines.append("  ('-' = not received; 'Nx' = N candidate licences on that channel, unidentified;\n"
                 "   TV level is channel power over the band floor, S/N is from the pilot or the power)")
    return "\n".join(lines)


if __name__ == "__main__":
    sys.exit(main())


def main_cli() -> None:
    sys.exit(main())
