"""Write a self-contained Leaflet map of the likelihood surface, transmitters and estimate."""
from __future__ import annotations

import json

import numpy as np

TEMPLATE = """<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>FMLocate estimate</title>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/leaflet/1.9.4/leaflet.min.css">
<script src="https://cdnjs.cloudflare.com/ajax/libs/leaflet/1.9.4/leaflet.min.js"></script>
<style>
html,body{margin:0;height:100%;font:14px -apple-system,system-ui,sans-serif}
#map{position:absolute;inset:0}
#info{position:absolute;z-index:1000;top:10px;right:10px;max-width:320px;background:#fffe;
  padding:10px 12px;border-radius:8px;box-shadow:0 1px 6px #0004;line-height:1.4}
#info b{font-size:15px}
</style></head><body>
<div id="map"></div><div id="info"></div>
<script>
const D = __DATA__;
const map = L.map('map');
L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png',
  {maxZoom: 15, attribution: '&copy; OpenStreetMap contributors'}).addTo(map);
// Heatmap: posterior probability, square-root scaled so the tails stay visible.
const g = D.grid, H = g.lats.length, W = g.lons.length;
const cv = document.createElement('canvas'); cv.width = W; cv.height = H;
const ctx = cv.getContext('2d'), img = ctx.createImageData(W, H);
for (let i = 0; i < H; i++) for (let j = 0; j < W; j++) {
  const v = g.v[i * W + j] / 255, k = ((H - 1 - i) * W + j) * 4;
  img.data[k] = 255; img.data[k + 1] = Math.round(220 * (1 - v)); img.data[k + 2] = 0;
  img.data[k + 3] = Math.round(200 * Math.sqrt(v));
}
ctx.putImageData(img, 0, 0);
const dy = (g.lats[1] - g.lats[0]) / 2, dx = (g.lons[1] - g.lons[0]) / 2;
const bounds = [[g.lats[0] - dy, g.lons[0] - dx], [g.lats[H - 1] + dy, g.lons[W - 1] + dx]];
L.imageOverlay(cv.toDataURL(), bounds, {opacity: 0.85}).addTo(map);
for (const t of D.tx) {
  L.circleMarker([t.lat, t.lon], {radius: 4 + Math.log10(1 + t.erp_kw) * 3, color: t.rx ? '#1565c0' : '#888',
    weight: 2, fillOpacity: t.rx ? 0.6 : 0.2})
   .bindTooltip(`${t.call} ${t.freq.toFixed(1)} (${t.service}) ${t.erp_kw} kW, ${t.haat_m ?? '?'} m HAAT` +
                (t.rx ? `<br>level ${t.level} dB` : '')).addTo(map);
}
const est = L.marker([D.est.lat, D.est.lon]).addTo(map).bindTooltip('Estimate (MAP)');
L.circle([D.est.lat, D.est.lon], {radius: D.r90 * 1000, color: '#c62828', fill: false, dashArray: '6 4'}).addTo(map);
if (D.truth) L.circleMarker([D.truth.lat, D.truth.lon], {radius: 7, color: '#2e7d32', weight: 3, fill: false})
  .bindTooltip('Known location').addTo(map);
map.fitBounds(L.latLngBounds([[D.est.lat, D.est.lon], ...D.tx.filter(t => t.rx).map(t => [t.lat, t.lon])]).pad(0.3));
document.getElementById('info').innerHTML = `<b>FMLocate</b><br>
Estimate ${D.est.lat.toFixed(3)}, ${D.est.lon.toFixed(3)}<br>
90% region &asymp; ${D.r90.toFixed(0)} km radius (dashed)<br>
${D.truth ? `Known location: ${D.truth.err.toFixed(1)} km away (green)<br>` : ''}
Blue: transmitters on received channels. Grey: other licences.<br>
<small>Heat = posterior probability of position. ${D.when}</small>`;
</script></body></html>
"""


def write(res: dict, meas: dict, path: str, truth: dict | None = None) -> None:
    g = res["grid"]
    v = np.round(255 * g["post"] / g["post"].max()).astype(int).ravel().tolist()
    levels = {c["freq"]: c["level_db"] for c in res["channels"] if c["detected"]}
    tx = []
    seen = set()
    for c in res["channels"]:
        if c["station"]:
            s = c["station"]
            seen.add((s["call"], c["freq"]))
            tx.append({**s, "freq": c["freq"], "service": "", "rx": c["detected"],
                       "level": levels.get(c["freq"])})
    for t in res["transmitters"]:
        if (t["call"], round(t["freq"], 1)) not in seen:
            tx.append({**t, "rx": False})
    data = {"grid": {"lats": g["lats"], "lons": g["lons"], "v": v}, "tx": tx, "est": res["estimate"],
            "r90": res["radius90_km"], "truth": truth, "when": meas.get("time_utc", "")}
    with open(path, "w") as fh:
        fh.write(TEMPLATE.replace("__DATA__", json.dumps(data)))
