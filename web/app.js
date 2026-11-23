// Page controller: talks to the Pyodide worker and draws the result.
// (Plotting/controls live here in JS; the physics lives in Python under sim/.)

const runBtn = document.getElementById("run");
const statusEl = document.getElementById("status");
const diagEl = document.getElementById("diagnostics");
const canvas = document.getElementById("map");
const ctx = canvas.getContext("2d");

const worker = new Worker("worker.js");

runBtn.addEventListener("click", () => {
  statusEl.textContent = "running…";
  runBtn.disabled = true;
  worker.postMessage({});
});

worker.onmessage = (e) => {
  const msg = e.data;
  if (!msg.ok) {
    statusEl.textContent = "error";
    diagEl.textContent = "Error: " + msg.error;
    return;
  }
  if (msg.ready) {
    // Pyodide + sim loaded: enable the button and run once automatically.
    statusEl.textContent = "ready";
    runBtn.disabled = false;
    worker.postMessage({});
    statusEl.textContent = "running…";
    return;
  }
  render(msg.data);
  statusEl.textContent = "done";
  runBtn.disabled = false;
};

// Simple blue -> green -> yellow -> red colour ramp for t in [0, 1].
function colormap(t) {
  t = Math.max(0, Math.min(1, t));
  const stops = [
    [0.0, [8, 29, 88]],
    [0.25, [34, 94, 168]],
    [0.5, [65, 182, 196]],
    [0.75, [254, 204, 92]],
    [1.0, [227, 26, 28]],
  ];
  for (let i = 1; i < stops.length; i++) {
    if (t <= stops[i][0]) {
      const [t0, c0] = stops[i - 1];
      const [t1, c1] = stops[i];
      const f = (t - t0) / (t1 - t0);
      return c0.map((c, k) => Math.round(c + f * (c1[k] - c)));
    }
  }
  return stops[stops.length - 1][1];
}

function render(data) {
  const { nx, nz, C, peak, Pe, Cr, interfaces, depth } = data;

  diagEl.textContent =
    `cell Peclet  Pe = ${Pe.toFixed(3)}` +
    (Pe > 2 ? "   ** WARNING: Pe > 2, central differences may oscillate **" : "") +
    `\ncell Courant Cr = ${Cr.toFixed(3)}` +
    `\npeak concentration = ${peak.toFixed(2)} mg/L`;

  // Draw the field at native grid resolution; CSS scales it up (pixelated).
  canvas.width = nx;
  canvas.height = nz;
  const img = ctx.createImageData(nx, nz);
  const scale = peak > 0 ? 1 / peak : 1;
  for (let ix = 0; ix < nx; ix++) {
    for (let iz = 0; iz < nz; iz++) {
      const v = C[ix * nz + iz];          // flattened row-major [ix, iz]
      const [r, g, b] = colormap(v * scale);
      const p = (iz * nx + ix) * 4;       // pixel (x=ix, y=iz)
      img.data[p] = r;
      img.data[p + 1] = g;
      img.data[p + 2] = b;
      img.data[p + 3] = 255;
    }
  }
  ctx.putImageData(img, 0, 0);

  // Dashed layer-interface lines at their depths.
  ctx.strokeStyle = "rgba(255,255,255,0.8)";
  ctx.lineWidth = 0.5;
  ctx.setLineDash([2, 2]);
  for (const z of interfaces) {
    const y = (z / depth) * nz;
    ctx.beginPath();
    ctx.moveTo(0, y);
    ctx.lineTo(nx, y);
    ctx.stroke();
  }
}
