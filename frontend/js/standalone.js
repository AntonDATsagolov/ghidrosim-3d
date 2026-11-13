/**
 * Интерфейс автономной HTML-версии.
 *
 * Работает без сервера и без Python: вычисления выполняет движок
 * engine.js, данные проекта хранятся в памяти (и выгружаются в JSON).
 * Разметка страницы общая с серверной версией (frontend/index.html).
 */
import {
  makeGrid, buildPropertyModels, buildSources, solveDarcy, traceStreamlines,
  seedsFromWells, stopCellsFromWells, simulate3D, calculateReserves,
  simulateFlowline, autoMapColumns, parseDelimited, rowsToWells,
} from './engine.js';
import { buildDemoProject } from './demo.js';
import { Viewer3D, cssColor } from './viewer3d.js';
import { lineChart, heatmap, drawLegend } from './charts.js';

const $ = (s) => document.querySelector(s);
const $$ = (s) => Array.from(document.querySelectorAll(s));

const state = {
  project: null,
  results: null,
  viewer: null,
  sliceCache: {},
  streamlines: [],
  layerTouched: false,
  editingWell: null,
  importRows: [],
  importMapping: {},
};

function toast(msg, kind = 'info') {
  const el = $('#toast');
  el.textContent = msg;
  el.className = `toast ${kind}`;
  clearTimeout(toast._t);
  toast._t = setTimeout(() => el.classList.add('hidden'), 4200);
}

function status(msg) { $('#statusBar').textContent = msg; }

function fmtNum(v) {
  const a = Math.abs(v);
  if (a >= 1000) return v.toFixed(0);
  if (a >= 10) return v.toFixed(1);
  if (a >= 1) return v.toFixed(2);
  if (a === 0) return '0';
  return v.toFixed(3);
}

// ------------------------------------------------------------------ вкладки

function initTabs() {
  $('#leftTabs').addEventListener('click', (e) => {
    const btn = e.target.closest('.tab');
    if (!btn) return;
    $$('#leftTabs .tab').forEach((b) => b.classList.remove('active'));
    btn.classList.add('active');
    $$('.panel-left .tabpane').forEach((p) => p.classList.remove('active'));
    $(`#tab-${btn.dataset.tab}`).classList.add('active');
  });
  $('#rightTabs').addEventListener('click', (e) => {
    const btn = e.target.closest('.tab');
    if (!btn) return;
    $$('#rightTabs .tab').forEach((b) => b.classList.remove('active'));
    btn.classList.add('active');
    $$('.panel-right .tabpane').forEach((p) => p.classList.remove('active'));
    $(`#rtab-${btn.dataset.rtab}`).classList.add('active');
    if (btn.dataset.rtab === 'slice') drawSliceMap();
    if (btn.dataset.rtab === 'flowline') selectFlowline();
  });
}

// ------------------------------------------------------------------ проект

function newProject() {
  state.project = buildDemoProject();
  state.results = null;
  state.sliceCache = {};
  state.streamlines = [];
  state.layerTouched = false;
  refresh();
  toast('Создан демонстрационный проект. Нажмите «Расчёт».');
}

function refresh() {
  const p = state.project;
  $('#projectLabel').textContent = `${p.name} · сетка ${p.grid.nx}×${p.grid.ny}×${p.grid.nz}`;
  renderWells();
  renderGridForm();
  renderChemForm();
  renderBCForm();
  state.viewer.setBlock(p.grid);
  state.viewer.setWells(p.wells, p.grid);
  updateLayerSlider();
  renderSummary();
}

function renderWells() {
  const tb = $('#wellTable tbody');
  tb.innerHTML = '';
  for (const w of state.project.wells) {
    const tr = document.createElement('tr');
    const scr = (w.screens && w.screens.length)
      ? `${w.screens[0].top}–${w.screens[0].bottom}` : '—';
    tr.innerHTML = `
      <td title="${w.name}">${w.name}</td>
      <td><span class="tag ${w.well_type}">${w.well_type === 'injection' ? 'зак' : w.well_type === 'production' ? 'отк' : 'набл'}</span></td>
      <td class="mono">${Math.round(w.trajectory[0].x)}</td>
      <td class="mono">${Math.round(w.trajectory[0].y)}</td>
      <td class="mono">${Math.round(w.trajectory[0].z)}</td>
      <td class="mono">${scr}</td>
      <td class="mono">${Number(w.rate).toFixed(1)}</td>`;
    tr.addEventListener('click', () => editWell(w));
    tb.appendChild(tr);
  }
  $('#wellCount').textContent = `${state.project.wells.length} скважин`;
}

function editWell(w) {
  state.editingWell = w;
  const f = $('#wellForm');
  f.classList.remove('hidden');
  $('#wellFormTitle').textContent = `Скважина: ${w.name}`;
  f.name.value = w.name;
  f.well_type.value = w.well_type;
  f.x.value = w.trajectory[0].x;
  f.y.value = w.trajectory[0].y;
  f.z.value = w.trajectory[0].z;
  f.acid_in.value = w.acid_in ?? 20;
  f.screen_top.value = w.screens?.[0]?.top ?? '';
  f.screen_bottom.value = w.screens?.[0]?.bottom ?? '';
  f.rate.value = w.rate ?? 0;
}

function renderGridForm() {
  const g = state.project.grid;
  const f = $('#gridForm');
  ['x0', 'y0', 'z0', 'dx', 'dy', 'dz', 'nx', 'ny', 'nz'].forEach((k) => { f[k].value = g[k]; });
  $('#cellCount').textContent = (g.nx * g.ny * g.nz).toLocaleString('ru');
}

function renderChemForm() {
  const c = state.project.chem;
  const f = $('#chemForm');
  ['acid_concentration', 'fe3_conc', 'feoh3_conc', 'caco3_conc', 'integration_days',
    'temperature', 'activation_energy', 'beta', 'k1', 'k2', 'k3', 'k4',
    'flowline_k1', 'flowline_k2', 'flowline_k3', 'flowline_k4', 'flowline_beta',
    'flowline_acid'].forEach((k) => { if (f[k] && c[k] != null) f[k].value = c[k]; });
}

const BC_FACES = [
  ['x_min', 'Грань Xmin (запад)'], ['x_max', 'Грань Xmax (восток)'],
  ['y_min', 'Грань Ymin (юг)'], ['y_max', 'Грань Ymax (север)'],
  ['z_min', 'Крыша блока'], ['z_max', 'Подошва блока'],
];

function renderBCForm() {
  const p = state.project;
  const rows = $('#bcRows');
  rows.innerHTML = '';
  for (const [face, label] of BC_FACES) {
    const cur = p.bcFaces?.[face] || 'neumann';
    const div = document.createElement('div');
    div.className = 'bc-row';
    div.innerHTML = `
      <span>${label}</span>
      <select data-face="${face}">
        <option value="neumann"${cur === 'neumann' ? ' selected' : ''}>Неймана</option>
        <option value="dirichlet"${cur === 'dirichlet' ? ' selected' : ''}>Дирихле</option>
      </select>
      <input type="number" step="0.5" placeholder="напор, м" data-head="${face}"
             value="${p.bcHeads?.[face] ?? 0}"${cur === 'dirichlet' ? '' : ' disabled'}>`;
    const sel = div.querySelector('select');
    const inp = div.querySelector('input');
    sel.addEventListener('change', () => { inp.disabled = sel.value !== 'dirichlet'; });
    rows.appendChild(div);
  }
}

// ------------------------------------------------------------------ расчёт

function runCalc() {
  const overlay = $('#progressOverlay');
  overlay.classList.remove('hidden');
  $('#progressText').textContent = 'Расчёт гидродинамики и массопереноса…';
  $('#progressBar').value = 10;
  status('Идёт расчёт…');

  setTimeout(() => {
    try {
      const t0 = performance.now();
      const res = computeAll();
      res.elapsed_sec = (performance.now() - t0) / 1000;
      state.results = res;
      state.sliceCache = {};
      refresh();
      renderCurves();
      renderSummary();
      loadSlice();
      loadStreamlines();
      $('#progressBar').value = 100;
      toast(`Расчёт выполнен за ${res.elapsed_sec.toFixed(1)} с`);
      status('Расчёт завершён');
    } catch (e) {
      console.error(e);
      toast(`Ошибка расчёта: ${e.message}`, 'error');
      status('Ошибка расчёта');
    } finally {
      overlay.classList.add('hidden');
    }
  }, 60);
}

function computeAll() {
  const p = state.project;
  const geo = readGeoForm();
  const grid = makeGrid(p.grid);
  const props = buildPropertyModels(p, geo);

  const { q: sources } = buildSources(grid, p.wells);
  const hydro = solveDarcy(grid, props.k, sources, props.porosity,
    p.bcFaces ?? {}, p.bcHeads ?? {}, geo.kv_ratio);

  const { pts, labels } = seedsFromWells(p.wells, grid);
  const stops = stopCellsFromWells(p.wells, grid);
  const lines = traceStreamlines(grid, hydro.vx, hydro.vy, hydro.vz, pts, labels, stops, {
    stagnationThreshold: geo.stagnation_threshold,
    maxResidence: Math.max(365, (p.chem.integration_days ?? 365) * 3),
  });

  const every = geo.output_every;
  const times = [];
  for (let t = every; t <= (p.chem.integration_days ?? 365) + 1e-9; t += every) times.push(t);
  const mt = simulate3D(grid, props.k, props.porosity, props.uranium, props.clay,
    hydro.faceX, hydro.faceY, hydro.faceZ, p.wells, p.chem, times);

  const reserves = calculateReserves(grid, mt.residual[mt.residual.length - 1],
    mt.acid[mt.acid.length - 1], { acidThreshold: geo.acid_threshold });

  return {
    grid, props, hydro, lines, mt, reserves,
    property_stats: {
      k_min: Math.min(...props.k), k_max: Math.max(...props.k),
      k_mean: props.k.reduce((a, b) => a + b, 0) / props.k.length,
      u_mean: props.uranium.reduce((a, b) => a + b, 0) / props.uranium.length,
    },
    hydro_summary: {
      residual: hydro.residual, total_in: hydro.totalIn, total_out: hydro.totalOut,
      boundary_in: hydro.boundaryIn, boundary_out: hydro.boundaryOut,
      head_min: Math.min(...hydro.head), head_max: Math.max(...hydro.head),
      speed_max: Math.max(...hydro.speed),
      balance_ok: Math.abs(hydro.residual) < Math.max(1e-6, 0.01 * Math.max(hydro.totalIn, hydro.totalOut, 1)),
    },
    curves: {
      time: mt.times, pr_uranium: mt.prUranium, vr_acid: mt.vrAcid,
      cumulative_u: mt.cumulativeU, cumulative_acid: mt.cumulativeAcid,
      prod_volume: mt.prodVolume, inj_volume: mt.injVolume,
    },
    reserves,
    stagnant_fraction: hydro.speed.reduce((a, b) => a + (b < geo.stagnation_threshold ? 1 : 0), 0) / hydro.speed.length,
  };
}

function readGeoForm() {
  const g = $('#geoForm');
  return {
    interpolation: g.interpolation.value,
    variogram_model: g.variogram_model.value,
    variogram_sill: parseFloat(g.variogram_sill.value) || 1,
    variogram_range: parseFloat(g.variogram_range.value) || 150,
    variogram_nugget: parseFloat(g.variogram_nugget.value) || 0,
    anisotropy: [
      parseFloat(g.aniso_x.value) || 1,
      parseFloat(g.aniso_y.value) || 1,
      parseFloat(g.aniso_z.value) || 1,
    ],
    kv_ratio: parseFloat(g.kv_ratio.value) || 0.2,
    output_every: Math.max(1, parseFloat($('#chemForm').integration_days.value) / 12 || 30),
    acid_threshold: 5,
    stagnation_threshold: 0.05,
  };
}

// ------------------------------------------------------------------ поля

const FIELD_META = {
  head: 'Напор, м', speed: 'Скорость, м/сут', k: 'Проницаемость, м/сут',
  porosity: 'Пористость', uranium: 'Уран в руде, кг/т',
  acid_final: 'Кислотность раствора, г/л', uranium_final: 'Уран в растворе, мг/л',
  residual_final: 'Остаточный уран, кг/т',
};

function fieldArray(name) {
  const r = state.results;
  if (!r) return null;
  switch (name) {
    case 'head': return r.hydro.head;
    case 'speed': return r.hydro.speed;
    case 'k': return r.props.k;
    case 'porosity': return r.props.porosity;
    case 'uranium': return r.props.uranium;
    case 'acid_final': return r.mt.acid[r.mt.acid.length - 1];
    case 'uranium_final': return r.mt.uranium[r.mt.uranium.length - 1];
    case 'residual_final': return r.mt.residual[r.mt.residual.length - 1];
    default: return null;
  }
}

function loadSlice() {
  const name = $('#fieldSelect').value;
  const arr = fieldArray(name);
  if (!arr) return;
  const grid = state.results.grid;
  const layer = parseInt($('#layerSlider').value, 10);
  const n = grid.n ?? grid.nx * grid.ny * grid.nz;
  const data = new Float64Array(n);
  let mn = Infinity, mx = -Infinity;
  for (let k = 0; k < grid.nz; k++) for (let j = 0; j < grid.ny; j++) for (let i = 0; i < grid.nx; i++) {
    const v = arr[i + grid.nx * (j + grid.ny * k)];
    if (k === layer) data[j * grid.nx + i] = v;
    if (v < mn) mn = v;
    if (v > mx) mx = v;
  }
  state.viewer.setSlice(data, grid.nx, grid.ny, grid, layer, {
    opacity: parseInt($('#sliceOpacity').value, 10) / 100, min: mn, max: mx,
  });
  $('#legMin').textContent = fmtNum(mn);
  $('#legMid').textContent = fmtNum((mn + mx) / 2);
  $('#legMax').textContent = fmtNum(mx);
  $('#legTitle').textContent = `${FIELD_META[name] || name} · слой ${layer + 1}`;
  drawLegend($('#legendCanvas'), mn, mx, (t) => cssColor(t));
  state.sliceStats = { min: mn, max: mx, data, layer };
  if ($('#rtab-slice').classList.contains('active')) drawSliceMap();
}

function updateLayerSlider() {
  const grid = state.project.grid;
  const sl = $('#layerSlider');
  sl.max = String(grid.nz - 1);
  if (!state.layerTouched) sl.value = String(Math.max(0, Math.floor(grid.nz / 2) - 1));
  sl.value = String(Math.min(parseInt(sl.value, 10) || 0, grid.nz - 1));
  $('#layerLabel').textContent = `${parseInt(sl.value, 10) + 1} / ${grid.nz}`;
}

function drawSliceMap() {
  if (!state.sliceStats) return;
  const g = state.project.grid;
  const { min, max, data, layer } = state.sliceStats;
  heatmap($('#sliceCanvas'), data, g.nx, g.ny, {
    norm: (v) => (max > min ? (v - min) / (max - min) : 0.5),
    color: (t) => cssColor(t),
    xMin: `${Math.round(g.x0)} м`, xMax: `${Math.round(g.x0 + g.nx * g.dx)} м`,
    yMin: `${Math.round(g.y0)} м`, yMax: `${Math.round(g.y0 + g.ny * g.dy)} м`,
    title: `Слой ${layer + 1} (z = ${Math.round(g.z0 + (layer + 0.5) * g.dz)} м)`,
  });
  $('#sliceInfo').textContent = `диапазон ${fmtNum(min)} … ${fmtNum(max)}`;
}

// ------------------------------------------------------------------ кривые

function renderCurves() {
  const c = state.results?.curves;
  if (!c) return;
  lineChart($('#curveChart'), [
    { values: c.pr_uranium, color: '#ff5252', x: c.time, axis: 'left' },
    { values: c.vr_acid, color: '#4d8dff', x: c.time, axis: 'right' },
  ], { dualAxis: true, xLabel: 'время, сут' });
  lineChart($('#extractChart'), [
    { values: c.cumulative_u, color: '#8bc34a', x: c.time, axis: 'left' },
  ], { xLabel: 'время, сут' });
}

function renderSummary() {
  const r = state.results;
  const box = $('#summaryBox');
  if (!r) {
    box.innerHTML = '<p class="muted">Нажмите «Расчёт», чтобы увидеть результаты.</p>';
    return;
  }
  const h = r.hydro_summary;
  const res = r.reserves;
  const ps = r.property_stats;
  const rows = [
    ['Закачка, м³/сут', fmtNum(h.total_in)],
    ['Отбор, м³/сут', fmtNum(h.total_out)],
    ['Приток через границы, м³/сут', fmtNum(h.boundary_in)],
    ['Отток через границы, м³/сут', fmtNum(h.boundary_out)],
    ['Невязка баланса, м³/сут', h.residual.toExponential(2)],
    ['Напор, м', `${fmtNum(h.head_min)} … ${fmtNum(h.head_max)}`],
    ['Макс. скорость, м/сут', fmtNum(h.speed_max)],
    ['Застойные зоны, %', ((r.stagnant_fraction || 0) * 100).toFixed(1)],
    ['Проницаемость, м/сут', `${fmtNum(ps.k_min)} … ${fmtNum(ps.k_max)}`],
    ['Среднее содержание U, кг/т', fmtNum(ps.u_mean)],
  ];
  const rrows = [
    ['Запас U в закисляемом объёме, т', fmtNum(res.reserveT)],
    ['Закисляемый объём, м³', fmtNum(res.acidifiedVolume)],
    ['Площадь, га', fmtNum(res.areaHa)],
    ['Средневзвешенное содержание, кг/т', fmtNum(res.avgGrade)],
    ['Средняя мощность, м', fmtNum(res.avgThickness)],
    ['Средневзвешенный метропроцент', fmtNum(res.weightedMetrePercent)],
  ];
  box.innerHTML = `
    <h3>Гидродинамика</h3>
    <table class="kv">${rows.map(([k, v]) => `<tr><td>${k}</td><td class="mono">${v}</td></tr>`).join('')}</table>
    <h3>Запасы</h3><table class="kv">${rrows.map(([k, v]) => `<tr><td>${k}</td><td class="mono">${v}</td></tr>`).join('')}</table>
    ${h.balance_ok ? '' : '<p class="warn">⚠ Невязка баланса превышает 1% — проверьте расходы скважин и граничные условия.</p>'}`;
}

// ------------------------------------------------------------------ линии тока

function loadStreamlines() {
  const lines = state.results?.lines ?? [];
  state.streamlines = lines;
  state.viewer.setStreamlines(lines, state.results.grid);
  const sel = $('#flowlineSelect');
  sel.innerHTML = '';
  lines.slice(0, 300).forEach((l, i) => {
    const o = document.createElement('option');
    o.value = String(i);
    o.textContent = `${l.wellIn || '—'} → ${l.wellOut || 'выход'} · ${Math.round(l.residence)} сут`;
    sel.appendChild(o);
  });
}

function selectFlowline() {
  const idx = parseInt($('#flowlineSelect').value || '0', 10);
  const lines = state.streamlines;
  if (!lines.length || !state.results) return;
  const l = lines[Math.min(idx, lines.length - 1)];
  const g = state.results.grid;
  const grade = state.results.props.uranium;
  const poro = state.results.props.porosity;
  const gradeVals = l.points.map((pt) => {
    const c = gridCentersSafe(g);
    const fx = clampIdx((pt[0] - g.x0) / g.dx - 0.5, g.nx);
    const fy = clampIdx((pt[1] - g.y0) / g.dy - 0.5, g.ny);
    const fz = clampIdx((pt[2] - g.z0) / g.dz - 0.5, g.nz);
    return sampleField(grade, g, fx, fy, fz);
  });
  const poroVals = l.points.map((pt) => {
    const fx = clampIdx((pt[0] - g.x0) / g.dx - 0.5, g.nx);
    const fy = clampIdx((pt[1] - g.y0) / g.dy - 0.5, g.ny);
    const fz = clampIdx((pt[2] - g.z0) / g.dz - 0.5, g.nz);
    return sampleField(poro, g, fx, fy, fz);
  });
  const prof = simulateFlowline(l.points, l.times, gradeVals, poroVals, state.project.chem);
  lineChart($('#flowlineChart'), [
    { values: prof.uranium, color: '#ff5252', x: prof.times, axis: 'left' },
    { values: prof.acid, color: '#4d8dff', x: prof.times, axis: 'left' },
    { values: prof.rel.map((v) => v * 100), color: '#8bc34a', x: prof.times, axis: 'left' },
  ], { xLabel: 'время пребывания, сут' });
  const wname = (id) => {
    const w = state.project.wells.find((x) => x.id === id);
    return w ? w.name : '';
  };
  $('#flowlineInfo').innerHTML = `
    <div>Линия: ${wname(l.wellIn) || '—'} → ${wname(l.wellOut) || 'выход из блока'}</div>
    <div>Время пребывания: ${fmtNum(l.residence)} сут</div>
    <div>Кислота на выходе: ${fmtNum(prof.acid[prof.acid.length - 1])} г/л</div>
    <div>Уран в растворе: ${fmtNum(prof.uranium[prof.uranium.length - 1])} мг/л</div>
    <div>Остаток U в руде: ${(100 * prof.rel[prof.rel.length - 1]).toFixed(1)} %</div>`;
}

function clampIdx(v, n) { return v < 0 ? 0 : (v > n - 1.000001 ? n - 1.000001 : v); }

function sampleField(field, g, fx, fy, fz) {
  const i0 = Math.floor(fx), j0 = Math.floor(fy), k0 = Math.floor(fz);
  const i1 = Math.min(i0 + 1, g.nx - 1), j1 = Math.min(j0 + 1, g.ny - 1), k1 = Math.min(k0 + 1, g.nz - 1);
  const tx = fx - i0, ty = fy - j0, tz = fz - k0;
  const at = (i, j, k) => field[i * g.ny * g.nz + j * g.nz + k];
  const c00 = at(i0, j0, k0) * (1 - tx) + at(i1, j0, k0) * tx;
  const c10 = at(i0, j1, k0) * (1 - tx) + at(i1, j1, k0) * tx;
  const c01 = at(i0, j0, k1) * (1 - tx) + at(i1, j0, k1) * tx;
  const c11 = at(i0, j1, k1) * (1 - tx) + at(i1, j1, k1) * tx;
  const c0 = c00 * (1 - ty) + c10 * ty;
  const c1 = c01 * (1 - ty) + c11 * ty;
  return c0 * (1 - tz) + c1 * tz;
}

function gridCentersSafe(g) { return null; }

// ------------------------------------------------------------------ импорт

function initImport() {
  $('#btnImport').addEventListener('click', () => $('#importModal').classList.remove('hidden'));
  $('#btnCloseImport').addEventListener('click', () => $('#importModal').classList.add('hidden'));
  $('#importFile').addEventListener('change', (e) => {
    const file = e.target.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      const text = String(reader.result);
      const { columns, data } = parseDelimited(text);
      if (!columns.length) { toast('Не удалось разобрать файл', 'error'); return; }
      state.importMapping = autoMapColumns(columns);
      state.importRows = rowsToWells(data, state.importMapping);
      renderMapping(columns);
      $('#importCount').textContent = `Распознано скважин: ${state.importRows.length}`;
      status('Файл загружен');
    };
    reader.readAsText(file, 'utf-8');
  });
  $('#btnApplyImport').addEventListener('click', () => {
    if (!state.importRows.length) { toast('Нет данных для загрузки', 'warn'); return; }
    for (const r of state.importRows) {
      const id = 'W' + Math.random().toString(36).slice(2, 8);
      const screens = [];
      if (r.screen_top != null && r.screen_bottom != null && r.screen_bottom > r.screen_top) {
        screens.push({ top: r.screen_top, bottom: r.screen_bottom });
      }
      state.project.wells.push({
        id, name: r.name || id,
        well_type: r.well_type || 'injection',
        trajectory: [{ x: r.x, y: r.y, z: r.z || 0 }],
        screens, screen: screens[0] ?? null,
        rate: r.rate || 0,
        acid_in: r.acid_in ?? 20,
      });
      const P = state.project.props;
      if (r.k != null) P.permeability.push([r.x, r.y, r.z || 0, r.k]);
      if (r.uranium != null) P.uranium.push([r.x, r.y, r.z || 0, r.uranium]);
      if (r.porosity != null) P.porosity.push([r.x, r.y, r.z || 0, r.porosity]);
      if (r.clay != null) P.clay.push([r.x, r.y, r.z || 0, r.clay]);
      if (r.lithology) P.lithology.push([r.x, r.y, r.z || 0, r.lithology]);
    }
    $('#importModal').classList.add('hidden');
    refresh();
    toast('Данные загружены в проект');
  });
}

function renderMapping(columns) {
  const box = $('#mappingRows');
  box.innerHTML = '';
  const labels = {
    name: 'Наименование', x: 'X', y: 'Y', z: 'Z (отметка)', k: 'Проницаемость',
    porosity: 'Пористость', uranium: 'Уран', clay: 'Глина', lithology: 'Литология',
    rate: 'Дебит', well_type: 'Тип', acid_in: 'Кислотность ВР',
    screen_top: 'Верх фильтра', screen_bottom: 'Низ фильтра',
  };
  for (const [field, label] of Object.entries(labels)) {
    const row = document.createElement('div');
    row.className = 'map-row';
    const opts = ['', ...columns].map((c) =>
      `<option value="${c}"${c === state.importMapping[field] ? ' selected' : ''}>${c || '—'}</option>`).join('');
    row.innerHTML = `<span>${label}</span><select data-field="${field}">${opts}</select>`;
    box.appendChild(row);
  }
  $('#importMapping').classList.remove('hidden');
  box.onchange = () => {
    box.querySelectorAll('select').forEach((s) => { state.importMapping[s.dataset.field] = s.value || null; });
    // пересчёт строк с новым сопоставлением
    const data = state.importRaw ?? [];
    state.importRows = rowsToWells(data, state.importMapping);
    $('#importCount').textContent = `Распознано скважин: ${state.importRows.length}`;
  };
}

// ------------------------------------------------------------------ сохранение

function download(filename, content, type = 'application/json') {
  const blob = new Blob([content], { type });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}

function saveProject() {
  download('project-ghidrosim.json', JSON.stringify(state.project, null, 1));
  toast('Проект сохранён в файл');
}

function openProject() {
  const input = document.createElement('input');
  input.type = 'file';
  input.accept = '.json';
  input.addEventListener('change', () => {
    const file = input.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      try {
        state.project = JSON.parse(String(reader.result));
        state.results = null;
        state.sliceCache = {};
        refresh();
        toast('Проект загружен');
      } catch (e) {
        toast('Не удалось прочитать файл', 'error');
      }
    };
    reader.readAsText(file, 'utf-8');
  });
  input.click();
}

function exportCsv() {
  const p = state.project;
  const lines = ['Скважина;X;Y;Z;Тип;Верх фильтра;Низ фильтра;Дебит'];
  for (const w of p.wells) {
    const s = w.screens?.[0];
    lines.push([w.name, w.trajectory[0].x, w.trajectory[0].y, w.trajectory[0].z,
      w.well_type, s ? s.top : '', s ? s.bottom : '', w.rate].join(';'));
  }
  download('skvazhiny.csv', '﻿' + lines.join('\n'), 'text/csv');
}

// ------------------------------------------------------------------ формы

function initForms() {
  $('#wellForm').addEventListener('submit', (e) => {
    e.preventDefault();
    const f = e.target;
    const item = {
      name: f.name.value, well_type: f.well_type.value,
      trajectory: [{ x: parseFloat(f.x.value), y: parseFloat(f.y.value), z: parseFloat(f.z.value) }],
      screens: [], screen: null,
      rate: parseFloat(f.rate.value), acid_in: parseFloat(f.acid_in.value),
    };
    if (f.screen_top.value !== '' && f.screen_bottom.value !== '') {
      item.screens = [{ top: parseFloat(f.screen_top.value), bottom: parseFloat(f.screen_bottom.value) }];
      item.screen = item.screens[0];
    }
    if (state.editingWell) {
      const i = state.project.wells.indexOf(state.editingWell);
      state.project.wells[i] = { ...state.editingWell, ...item };
    } else {
      item.id = 'W' + Math.random().toString(36).slice(2, 8);
      state.project.wells.push(item);
    }
    state.editingWell = null;
    f.classList.add('hidden');
    refresh();
    toast('Скважины сохранены');
  });

  $('#btnCancelWell').addEventListener('click', () => {
    state.editingWell = null;
    $('#wellForm').classList.add('hidden');
  });
  $('#btnAddWell').addEventListener('click', () => {
    state.editingWell = null;
    const f = $('#wellForm');
    f.reset();
    f.classList.remove('hidden');
    $('#wellFormTitle').textContent = 'Новая скважина';
  });

  $('#gridForm').addEventListener('submit', (e) => {
    e.preventDefault();
    const f = e.target;
    const g = state.project.grid;
    ['x0', 'y0', 'z0', 'dx', 'dy', 'dz', 'nx', 'ny', 'nz'].forEach((k) => {
      g[k] = parseFloat(f[k].value);
    });
    g.nx = Math.max(2, Math.round(g.nx));
    g.ny = Math.max(2, Math.round(g.ny));
    g.nz = Math.max(1, Math.round(g.nz));
    refresh();
    toast('Сетка обновлена');
  });

  $('#chemForm').addEventListener('submit', (e) => {
    e.preventDefault();
    const f = e.target;
    const c = state.project.chem;
    ['acid_concentration', 'fe3_conc', 'feoh3_conc', 'caco3_conc', 'integration_days',
      'temperature', 'activation_energy', 'beta', 'k1', 'k2', 'k3', 'k4',
      'flowline_k1', 'flowline_k2', 'flowline_k3', 'flowline_k4', 'flowline_beta',
      'flowline_acid'].forEach((k) => {
      if (f[k] && f[k].value !== '') c[k] = parseFloat(f[k].value);
    });
    refresh();
    toast('Параметры кинетики сохранены');
  });

  $('#bcForm').addEventListener('submit', (e) => {
    e.preventDefault();
    const p = state.project;
    p.bcFaces = {}; p.bcHeads = {};
    $$('#bcRows select').forEach((s) => { p.bcFaces[s.dataset.face] = s.value; });
    $$('#bcRows input').forEach((i) => { p.bcHeads[i.dataset.head] = parseFloat(i.value) || 0; });
    refresh();
    toast('Граничные условия применены');
  });

  $('#geoForm').addEventListener('submit', (e) => { e.preventDefault(); toast('Настройки сохранены'); });
}

// ------------------------------------------------------------------ init

function initViewerControls() {
  $('#fieldSelect').addEventListener('change', () => { if (state.results) loadSlice(); });
  $('#layerSlider').addEventListener('input', () => {
    state.layerTouched = true;
    $('#layerLabel').textContent = `${parseInt($('#layerSlider').value, 10) + 1} / ${state.project.grid.nz}`;
    if (state.results) loadSlice();
  });
  $('#sliceOpacity').addEventListener('input', () => {
    state.viewer.setSliceOpacity(parseInt($('#sliceOpacity').value, 10) / 100);
  });
  $('#showWells').addEventListener('change', (e) => state.viewer.setVisible('wells', e.target.checked));
  $('#showLines').addEventListener('change', (e) => state.viewer.setVisible('lines', e.target.checked));
  $('#showSlice').addEventListener('change', (e) => state.viewer.setVisible('slice', e.target.checked));
  $('#btnResetView').addEventListener('click', () => state.viewer.fitCamera());
  $('#flowlineSelect').addEventListener('change', selectFlowline);
  addEventListener('resize', () => state.viewer.resize());
}

function init() {
  initTabs();
  state.viewer = new Viewer3D($('#viewer'));
  initForms();
  initImport();
  initViewerControls();

  $('#btnDemo').addEventListener('click', newProject);
  $('#btnNew').addEventListener('click', newProject);
  $('#btnRun').addEventListener('click', runCalc);
  $('#btnSave').addEventListener('click', saveProject);
  $('#btnOpen').addEventListener('click', openProject);
  $('#btnExport').addEventListener('click', exportCsv);

  state.project = buildDemoProject();
  refresh();
  runCalc();
}

function start() {
  try {
    init();
  } catch (e) {
    console.error(e);
    const bar = document.querySelector('#statusBar');
    if (bar) bar.textContent = 'Ошибка инициализации: ' + (e && e.message ? e.message : e);
    const box = document.querySelector('#summaryBox');
    if (box) box.innerHTML = '<p class="warn">Ошибка: ' + (e && e.message ? e.message : e) + '</p>';
  }
}

if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', start);
} else {
  start();
}
