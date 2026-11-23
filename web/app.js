// Page controller: builds the input form, talks to the Pyodide worker, and
// draws the concentration-map animation and breakthrough curve.
// (All plotting/controls are plain JS here; the physics lives in Python/sim.)

// ----------------------------- DOM references -----------------------------
const $ = (id) => document.getElementById(id);
const statusEl = $("status");
const runBtn = $("run");
const diagEl = $("diagnostics");
const mapCanvas = $("map");
const mapCtx = mapCanvas.getContext("2d");
const cbarCanvas = $("colorbar");
const btCanvas = $("breakthrough");
const btCtx = btCanvas.getContext("2d");
const playBtn = $("play");
const frameSlider = $("frame");
const frameTimeEl = $("frame_time");

// ------------------------- Layer-table definitions ------------------------
// Columns shown for each velocity-input mode.
const COLS = {
  velocity: [
    { key: "thickness", label: "thickness [m]" },
    { key: "velocity", label: "velocity [m/day]" },
    { key: "dispersion", label: "dispersion [m²/day]" },
  ],
  darcy: [
    { key: "thickness", label: "thickness [m]" },
    { key: "K", label: "K [m/day]" },
    { key: "gradient", label: "gradient i [-]" },
    { key: "porosity", label: "porosity n [-]" },
    { key: "dispersion", label: "dispersion [m²/day]" },
  ],
};

// Separate stored rows per mode so switching modes does not lose data.
const layersByMode = {
  velocity: [
    { thickness: 10, velocity: 0.5, dispersion: 2.0 },
    { thickness: 10, velocity: 0.05, dispersion: 0.2 },
    { thickness: 20, velocity: 0.3, dispersion: 1.0 },
  ],
  darcy: [
    { thickness: 10, K: 5.0, gradient: 0.01, porosity: 0.1, dispersion: 2.0 },
    { thickness: 10, K: 0.5, gradient: 0.01, porosity: 0.1, dispersion: 0.2 },
    { thickness: 20, K: 3.0, gradient: 0.01, porosity: 0.1, dispersion: 1.0 },
  ],
};

let mode = "velocity";

function currentMode() {
  return document.querySelector('input[name="mode"]:checked').value;
}

// Read the row inputs back into the stored layers for a given mode.
function syncLayers(m) {
  const out = [];
  for (const rowEl of $("layer-rows").children) {
    const row = {};
    rowEl.querySelectorAll("input").forEach((inp) => { row[inp.dataset.key] = parseFloat(inp.value); });
    out.push(row);
  }
  layersByMode[m] = out;
}

// Depth is NOT an input: it is the sum of the layer thicknesses. Update the
// read-only display and return the total.
function updateDepth() {
  let d = 0;
  for (const inp of $("layer-rows").querySelectorAll('input[data-key="thickness"]')) {
    const v = parseFloat(inp.value);
    if (!Number.isNaN(v)) d += v;
  }
  $("depth_display").textContent = d > 0 ? (+d.toFixed(4)) + " m" : "—";
  return d;
}

function renderLayers() {
  const cols = COLS[mode];
  // Header labels aligned to the input grid, plus a spacer over the × column.
  const head = $("layer-head");
  head.innerHTML =
    `<div class="layer-grid" style="--ncols:${cols.length}">` +
    cols.map((c) => `<span>${c.label}</span>`).join("") +
    `</div><span class="del-spacer"></span>`;

  const rows = $("layer-rows");
  rows.innerHTML = "";
  layersByMode[mode].forEach((row, idx) => {
    const rowEl = document.createElement("div");
    rowEl.className = "layer-row";
    const grid = document.createElement("div");
    grid.className = "layer-grid";
    grid.style.setProperty("--ncols", cols.length);
    grid.innerHTML = cols.map((c) => {
      const v = row[c.key] !== undefined ? row[c.key] : "";
      return `<input type="number" step="any" data-key="${c.key}" value="${v}">`;
    }).join("");
    const del = document.createElement("button");
    del.type = "button";
    del.className = "del";
    del.title = "delete layer";
    del.textContent = "×";
    del.onclick = () => {
      syncLayers(mode);
      layersByMode[mode].splice(idx, 1);
      renderLayers();
    };
    rowEl.appendChild(grid);
    rowEl.appendChild(del);
    rows.appendChild(rowEl);
  });

  $("csv_hint").textContent =
    "CSV columns: " + cols.map((c) => c.key).join(", ") + " (one row per layer; a header row is ignored)";
  updateDepth();
}

// Recompute the depth readout whenever a thickness is edited.
$("layer-rows").addEventListener("input", updateDepth);

$("add_layer").onclick = () => {
  syncLayers(mode);
  const blank = {};
  COLS[mode].forEach((c) => (blank[c.key] = 0));
  layersByMode[mode].push(blank);
  renderLayers();
};

for (const r of document.querySelectorAll('input[name="mode"]')) {
  r.onchange = () => {
    syncLayers(mode);         // save edits from the old mode's rows
    mode = currentMode();
    renderLayers();
  };
}

// ------------------------------- CSV loading ------------------------------
$("csv").addEventListener("change", (e) => {
  const file = e.target.files[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = () => {
    const cols = COLS[mode];
    const rows = [];
    for (const line of reader.result.split(/\r?\n/)) {
      const cells = line.split(",").map((s) => s.trim());
      if (cells.length < cols.length) continue;
      const nums = cells.slice(0, cols.length).map(Number);
      if (nums.some((n) => Number.isNaN(n))) continue;   // skip header / bad rows
      const row = {};
      cols.forEach((c, i) => (row[c.key] = nums[i]));
      rows.push(row);
    }
    if (rows.length) {
      layersByMode[mode] = rows;
      renderLayers();
      statusEl.textContent = `loaded ${rows.length} layers from CSV`;
    } else {
      statusEl.textContent = "no valid rows found in CSV";
    }
    $("csv").value = "";   // allow re-loading the same file
  };
  reader.readAsText(file);
});

// Show/hide the leak-duration input based on the "continuous" checkbox.
$("continuous").onchange = () => {
  $("leak_label").style.display = $("continuous").checked ? "none" : "flex";
};
$("continuous").onchange();

// --------------------------- Gather parameters ----------------------------
function num(id) {
  return parseFloat($(id).value);
}

function gatherParams() {
  syncLayers(mode);
  const depth = layersByMode[mode].reduce((s, r) => s + (r.thickness || 0), 0);
  return {
    length: num("length"), depth: depth, dx: num("dx"), dz: num("dz"),
    dt: num("dt"), total_time: num("total_time"),
    velocity_mode: mode,
    layers: layersByMode[mode],
    source: {
      x: num("src_x"), width: num("src_width"), conc: num("src_conc"),
      leak_duration: $("continuous").checked ? null : num("leak_duration"),
    },
    receptor: { x: num("recv_x"), z: num("recv_z") },
    nframes: 40,
  };
}

// --------------------------- Stability pre-check --------------------------
// Crank–Nicolson ADI is unconditionally stable (no Δt limit), so the real
// failure mode is NOT blow-up but spurious oscillations when the cell Péclet
// number Pe = |v|·Δx/D exceeds 2. We compute the worst-case Pe and Cr from the
// inputs and warn before running.
function stabilityDiagnose(p) {
  let maxPe = 0, maxCr = 0;
  for (const row of p.layers) {
    const v = p.velocity_mode === "darcy"
      ? (row.K * row.gradient) / row.porosity
      : row.velocity;
    if (row.dispersion > 0) maxPe = Math.max(maxPe, Math.abs(v) * p.dx / row.dispersion);
    maxCr = Math.max(maxCr, Math.abs(v) * p.dt / p.dx);
  }
  return { maxPe, maxCr };
}

// ------------------------------- Worker I/O -------------------------------
const worker = new Worker("worker.js");

runBtn.onclick = () => {
  const params = gatherParams();
  const { maxPe, maxCr } = stabilityDiagnose(params);

  if (maxPe > 2) {
    const proceed = confirm(
      "⚠ Likely oscillatory (non-physical) solution.\n\n" +
      `Cell Péclet number Pe = ${maxPe.toFixed(2)}  (> 2)  at Δx = ${params.dx} m.\n\n` +
      "Crank–Nicolson is unconditionally stable, so it will NOT blow up, but the " +
      "central-difference advection term produces oscillations — wiggles and even " +
      "negative concentrations — whenever Pe > 2.\n" +
      "Fix: reduce Δx, or increase the dispersion D of the offending layer.\n\n" +
      `(Courant number Cr = ${maxCr.toFixed(2)} affects accuracy only, not stability.)\n\n` +
      "Run anyway?"
    );
    if (!proceed) { statusEl.textContent = "cancelled — reduce Δx or raise D (Pe > 2)"; return; }
  }

  stopAnimation();
  statusEl.textContent = "running…";
  runBtn.disabled = true;
  worker.postMessage(params);
};

worker.onmessage = (e) => {
  const msg = e.data;
  if (!msg.ok) {
    statusEl.textContent = "error";
    diagEl.textContent = "Error: " + msg.error;
    runBtn.disabled = false;
    return;
  }
  if (msg.ready) {
    statusEl.textContent = "ready";
    runBtn.disabled = false;
    return;
  }
  showResult(msg.data);
  statusEl.textContent = "done";
  runBtn.disabled = false;
};

// ------------------------- Results & animation ----------------------------
let R = null;          // latest result
let frameIndex = 0;
let timer = null;

function showResult(data) {
  R = data;
  frameIndex = 0;
  frameSlider.max = data.frames.length - 1;
  frameSlider.value = 0;
  frameSlider.disabled = false;
  playBtn.disabled = false;

  const warn = data.Pe > 2 ? '   <span class="warn">** Pe &gt; 2: central differences may oscillate **</span>' : "";
  diagEl.innerHTML =
    `cell Péclet  Pe = ${data.Pe.toFixed(3)}${warn}\n` +
    `cell Courant Cr = ${data.Cr.toFixed(3)}\n` +
    `receptor at x = ${data.receptor.x.toFixed(1)} m, z = ${data.receptor.z.toFixed(1)} m\n` +
    `peak concentration (colour scale max) = ${data.peak.toFixed(2)} mg/L`;

  drawColorbar(data.peak);
  positionMarkers(data);
  renderFrame(0);
}

frameSlider.oninput = () => {
  stopAnimation();
  renderFrame(parseInt(frameSlider.value, 10));
};

playBtn.onclick = () => {
  if (timer) { stopAnimation(); return; }
  playBtn.textContent = "❚❚ Pause";
  timer = setInterval(() => {
    frameIndex = (frameIndex + 1) % R.frames.length;
    frameSlider.value = frameIndex;
    renderFrame(frameIndex);
  }, 160);
};

function stopAnimation() {
  if (timer) { clearInterval(timer); timer = null; }
  playBtn.textContent = "▶ Play";
}

function renderFrame(i) {
  frameIndex = i;
  renderMap(R.frames[i]);
  frameTimeEl.textContent = `t = ${R.frame_times[i].toFixed(1)} days`;
  drawBreakthrough(R.frame_times[i]);
}

// --------------------------------- Drawing --------------------------------
// Viridis colormap (perceptually uniform), 11 samples interpolated linearly.
const VIRIDIS = [
  [68, 1, 84], [72, 40, 120], [62, 74, 137], [49, 104, 142], [38, 130, 142],
  [31, 158, 137], [53, 183, 121], [109, 205, 89], [180, 222, 44], [221, 228, 30], [253, 231, 37],
];
function colormap(t) {
  t = Math.max(0, Math.min(1, t));
  const x = t * (VIRIDIS.length - 1);
  const i = Math.min(Math.floor(x), VIRIDIS.length - 2);
  const f = x - i;
  const a = VIRIDIS[i], b = VIRIDIS[i + 1];
  return [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])].map(Math.round);
}

function renderMap(frame) {
  const { nx, nz, peak, interfaces, depth } = R;
  mapCanvas.width = nx;
  mapCanvas.height = nz;
  const img = mapCtx.createImageData(nx, nz);
  const scale = peak > 0 ? 1 / peak : 1;
  for (let ix = 0; ix < nx; ix++) {
    for (let iz = 0; iz < nz; iz++) {
      const [r, g, b] = colormap(frame[ix * nz + iz] * scale);  // flat [ix*nz+iz]
      const p = (iz * nx + ix) * 4;                              // pixel (x=ix, y=iz)
      img.data[p] = r; img.data[p + 1] = g; img.data[p + 2] = b; img.data[p + 3] = 255;
    }
  }
  mapCtx.putImageData(img, 0, 0);
  // Dashed layer interfaces at their depths.
  mapCtx.strokeStyle = "rgba(255,255,255,0.85)";
  mapCtx.lineWidth = 0.4;
  mapCtx.setLineDash([2, 2]);
  for (const z of interfaces) {
    const y = (z / depth) * nz;
    mapCtx.beginPath(); mapCtx.moveTo(0, y); mapCtx.lineTo(nx, y); mapCtx.stroke();
  }
}

function drawColorbar(peak) {
  const ctx = cbarCanvas.getContext("2d");
  const w = cbarCanvas.width, h = cbarCanvas.height;
  ctx.clearRect(0, 0, w, h);          // clear so repeated runs don't stack text

  const barH = 16;
  for (let x = 0; x < w; x++) {
    const [r, g, b] = colormap(x / (w - 1));
    ctx.fillStyle = `rgb(${r},${g},${b})`;
    ctx.fillRect(x, 0, 1, barH);
  }

  // Five evenly spaced ticks; end labels are left/right aligned so large
  // numbers (e.g. 1000) never overlap or run off the canvas.
  ctx.fillStyle = "#000";
  ctx.strokeStyle = "#000";
  ctx.font = "11px system-ui, sans-serif";
  const fmt = (v) => v === 0 ? "0"
    : (peak >= 100 ? Math.round(v).toString() : Number(v.toPrecision(3)).toString());
  const ticks = [0, 0.25, 0.5, 0.75, 1];
  ticks.forEach((fr, i) => {
    const x = fr * (w - 1);
    ctx.beginPath(); ctx.moveTo(x, barH); ctx.lineTo(x, barH + 3); ctx.stroke();
    ctx.textAlign = i === 0 ? "left" : (i === ticks.length - 1 ? "right" : "center");
    ctx.fillText(fmt(peak * fr), x, barH + 15);
  });
  ctx.textAlign = "center";
  ctx.fillText("concentration C [mg/L]", w / 2, h - 2);
}

// Overlay the source patch (bar) and receptor (dot) on the map, positioned by
// percentage so they track the responsive canvas size. Set once per result.
function positionMarkers(d) {
  const sm = $("src-marker"), rm = $("recv-marker");
  sm.style.left = ((d.source.x - d.source.width / 2) / d.length * 100) + "%";
  sm.style.width = (d.source.width / d.length * 100) + "%";
  sm.style.top = "0";
  sm.innerHTML = '<span class="mlabel">source</span>';
  rm.style.left = (d.receptor.x / d.length * 100) + "%";
  rm.style.top = (d.receptor.z / d.depth * 100) + "%";
  rm.innerHTML = '<span class="mlabel">receptor</span>';
}

function drawBreakthrough(markTime) {
  const times = R.times, vals = R.breakthrough;
  const W = btCanvas.width, H = btCanvas.height;
  const L = 55, Rm = 15, T = 15, B = 40;
  const pw = W - L - Rm, ph = H - T - B;
  btCtx.clearRect(0, 0, W, H);

  const xmax = times[times.length - 1] || 1;
  let ymax = 0;
  for (const v of vals) if (v !== null && v > ymax) ymax = v;
  ymax = ymax > 0 ? ymax * 1.1 : 1;

  const X = (t) => L + (t / xmax) * pw;
  const Y = (v) => T + ph - (v / ymax) * ph;

  // Axes
  btCtx.strokeStyle = "#333"; btCtx.lineWidth = 1; btCtx.setLineDash([]);
  btCtx.beginPath(); btCtx.moveTo(L, T); btCtx.lineTo(L, T + ph); btCtx.lineTo(L + pw, T + ph); btCtx.stroke();

  btCtx.fillStyle = "#333"; btCtx.font = "11px system-ui, sans-serif";
  btCtx.textAlign = "right";
  for (let k = 0; k <= 4; k++) {
    const v = (ymax / 4) * k, y = Y(v);
    btCtx.strokeStyle = "#eee"; btCtx.beginPath(); btCtx.moveTo(L, y); btCtx.lineTo(L + pw, y); btCtx.stroke();
    btCtx.fillText(v.toFixed(v < 10 ? 2 : 1), L - 5, y + 3);
  }
  btCtx.textAlign = "center";
  for (let k = 0; k <= 5; k++) {
    const t = (xmax / 5) * k;
    btCtx.fillText(t.toFixed(0), X(t), T + ph + 15);
  }
  btCtx.fillText("time [days]", L + pw / 2, H - 4);
  btCtx.save();
  btCtx.translate(12, T + ph / 2); btCtx.rotate(-Math.PI / 2);
  btCtx.fillText("C at receptor [mg/L]", 0, 0);
  btCtx.restore();

  // Curve
  btCtx.strokeStyle = "#1b6"; btCtx.lineWidth = 1.8; btCtx.beginPath();
  let started = false;
  for (let i = 0; i < times.length; i++) {
    if (vals[i] === null) continue;
    const x = X(times[i]), y = Y(vals[i]);
    if (!started) { btCtx.moveTo(x, y); started = true; } else btCtx.lineTo(x, y);
  }
  btCtx.stroke();

  // Marker at the currently displayed map time
  const mx = X(markTime);
  btCtx.strokeStyle = "#b00"; btCtx.setLineDash([4, 3]); btCtx.beginPath();
  btCtx.moveTo(mx, T); btCtx.lineTo(mx, T + ph); btCtx.stroke();
  btCtx.setLineDash([]);
}

// --------------------------------- Startup --------------------------------
renderLayers();
