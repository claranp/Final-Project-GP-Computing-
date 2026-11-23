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

// Read the table inputs back into the stored rows for a given mode.
function syncTableToState(m) {
  const rows = [];
  for (const tr of $("layers").tBodies[0].rows) {
    const row = {};
    COLS[m].forEach((col, i) => {
      row[col.key] = parseFloat(tr.cells[i].querySelector("input").value);
    });
    rows.push(row);
  }
  layersByMode[m] = rows;
}

function renderTable() {
  const cols = COLS[mode];
  const thead = $("layers").tHead;
  const tbody = $("layers").tBodies[0];
  thead.innerHTML = "<tr>" + cols.map((c) => `<th>${c.label}</th>`).join("") + "<th></th></tr>";
  tbody.innerHTML = "";
  for (const row of layersByMode[mode]) {
    const tr = tbody.insertRow();
    for (const col of cols) {
      const td = tr.insertCell();
      const val = row[col.key] !== undefined ? row[col.key] : "";
      td.innerHTML = `<input type="number" step="any" value="${val}">`;
    }
    const del = tr.insertCell();
    del.innerHTML = '<button type="button" class="del" title="remove">×</button>';
    del.querySelector("button").onclick = () => {
      syncTableToState(mode);
      layersByMode[mode].splice([...tbody.rows].indexOf(tr), 1);
      renderTable();
    };
  }
  $("csv_hint").textContent =
    "CSV columns: " + cols.map((c) => c.key).join(", ") + " (one row per layer; a header row is ignored)";
}

$("add_layer").onclick = () => {
  syncTableToState(mode);
  const blank = {};
  COLS[mode].forEach((c) => (blank[c.key] = 0));
  layersByMode[mode].push(blank);
  renderTable();
};

for (const r of document.querySelectorAll('input[name="mode"]')) {
  r.onchange = () => {
    syncTableToState(mode);   // save edits from the old mode's table
    mode = currentMode();
    renderTable();
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
      renderTable();
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
  syncTableToState(mode);
  return {
    length: num("length"), depth: num("depth"), dx: num("dx"), dz: num("dz"),
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

// ------------------------------- Worker I/O -------------------------------
const worker = new Worker("worker.js");

runBtn.onclick = () => {
  stopAnimation();
  statusEl.textContent = "running…";
  runBtn.disabled = true;
  worker.postMessage(gatherParams());
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
// Blue -> green -> yellow -> red colour ramp for t in [0, 1].
function colormap(t) {
  t = Math.max(0, Math.min(1, t));
  const stops = [
    [0.0, [8, 29, 88]], [0.25, [34, 94, 168]], [0.5, [65, 182, 196]],
    [0.75, [254, 204, 92]], [1.0, [227, 26, 28]],
  ];
  for (let i = 1; i < stops.length; i++) {
    if (t <= stops[i][0]) {
      const [t0, c0] = stops[i - 1], [t1, c1] = stops[i];
      const f = (t - t0) / (t1 - t0);
      return c0.map((c, k) => Math.round(c + f * (c1[k] - c)));
    }
  }
  return stops[stops.length - 1][1];
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
  const w = cbarCanvas.width, h = 18;
  for (let x = 0; x < w; x++) {
    const [r, g, b] = colormap(x / (w - 1));
    ctx.fillStyle = `rgb(${r},${g},${b})`;
    ctx.fillRect(x, 0, 1, h);
  }
  ctx.fillStyle = "#000";
  ctx.font = "11px system-ui, sans-serif";
  ctx.fillText("0", 0, h + 13);
  ctx.fillText((peak / 2).toFixed(1), w / 2 - 10, h + 13);
  ctx.fillText(peak.toFixed(1) + " mg/L", w - 60, h + 13);
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
renderTable();
