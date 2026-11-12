/** Главный модуль приложения: состояние, API, привязка интерфейса. */
import { Viewer3D, cssColor, colorFor } from './viewer3d.js';
import { lineChart, heatmap, drawLegend } from './charts.js';

const API = '';
const state = {
  pid: null,
  project: null,
  results: null,
  viewer: null,
  sliceCache: {},
  flowlines: [],
  editingWell: null,
};

const $ = (sel) => document.querySelector(sel);
const $$ = (sel) => Array.from(document.querySelectorAll(sel));

function toast(msg, kind = 'info') {
  const el = $('#toast');
  el.textContent = msg;
  el.className = `toast ${kind}`;
  clearTimeout(toast._t);
  toast._t = setTimeout(() => el.classList.add('hidden'), 4200);
}

function status(msg) {
  $('#statusBar').textContent = msg;
}

async function api(path, options = {}) {
  const res = await fetch(API + path, {
    headers: { 'Content-Type': 'application/json' },
    ...options,
  });
  if (!res.ok) {
    let detail = res.statusText;
    try { detail = (await res.json()).detail || detail; } catch (_) { /* noop */ }
    throw new Error(detail);
  }
  return res.json();
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

async function loadDemo() {
  status('Загрузка демо-проекта…');
  const p = await api('/api/demo');
  state.pid = p.id;
  await refreshProject();
  toast('Демо-блок загружен. Нажмите «Расчёт».');
}

async function newProject() {
  const p = await api('/api/projects', {
    method: 'POST',
    body: JSON.stringify({ name: 'Новый проект' }),
  });
  state.pid = p.id;
  await refreshProject();
  toast('Создан новый проект');
}

async function refreshProject() {
  const p = await api(`/api/projects/${state.pid}`);
  state.project = p;
  state.results = {
    hydro: p.hydro, reserves: p.reserves, curves: p.curves,
    property_stats: p.property_stats, stagnant_fraction: p.stagnant_fraction,
  };
  $('#projectLabel').textContent = `${p.name} · сетка ${p.grid.nx}×${p.grid.ny}×${p.grid.nz}`;
  renderWells();
  renderGridForm();
  renderChemForm();
  renderBCForm();
  state.viewer.setBlock(p.grid);
  state.viewer.setWells(p.wells, p.grid);
  state.sliceCache = {};
  updateLayerSlider();
  renderSummary();
  status('Готов к работе');
}

function renderWells() {
  const tb = $('#wellTable tbody');
  tb.innerHTML = '';
  for (const w of state.project.wells) {
    const tr = document.createElement('tr');
    const scr = w.screen_top != null ? `${w.screen_top}–${w.screen_bottom}` : '—';
    tr.innerHTML = `
      <td title="${w.name}">${w.name}</td>
      <td><span class="tag ${w.well_type}">${w.well_type === 'injection' ? 'зак' : w.well_type === 'production' ? 'отк' : 'набл'}</span></td>
      <td class="mono">${Math.round(w.x)}</td>
      <td class="mono">${Math.round(w.y)}</td>
      <td class="mono">${Math.round(w.z ?? 0)}</td>
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
  f.x.value = w.x; f.y.value = w.y; f.z.value = w.z ?? 0;
  f.acid_in.value = w.acid_in ?? 20;
  f.screen_top.value = w.screen_top ?? '';
  f.screen_bottom.value = w.screen_bottom ?? '';
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
    'flowline_acid'].forEach((k) => {
    if (f[k] && c[k] != null) f[k].value = c[k];
  });
}

const FACES = [
  ['x_min', 'Грань Xmin (запад)'], ['x_max', 'Грань Xmax (восток)'],
  ['y_min', 'Грань Ymin (юг)'], ['y_max', 'Грань Ymax (север)'],
  ['z_min', 'Крыша блока'], ['z_max', 'Подошва блока'],
];

function renderBCForm() {
  const p = state.project;
  const rows = $('#bcRows');
  rows.innerHTML = '';
  for (const [face, label] of FACES) {
    const cur = p.bc_faces?.[face] || p.boundary || 'neumann';
    const div = document.createElement('div');
    div.className = 'bc-row';
    div.innerHTML = `
      <span>${label}</span>
      <select data-face="${face}">
        <option value="neumann"${cur === 'neumann' ? ' selected' : ''}>Неймана</option>
        <option value="dirichlet"${cur === 'dirichlet' ? ' selected' : ''}>Дирихле</option>
      </select>
      <input type="number" step="0.5" placeholder="напор, м" data-head="${face}"
             value="${p.bc_heads?.[face] ?? p.boundary_head ?? 0}"${cur === 'dirichlet' ? '' : ' disabled'}>`;
    const sel = div.querySelector('select');
    const inp = div.querySelector('input');
    sel.addEventListener('change', () => { inp.disabled = sel.value !== 'dirichlet'; });
    rows.appendChild(div);
  }
}

// ------------------------------------------------------------------ расчёт

async function runCalc() {
  if (!state.pid) { toast('Сначала создайте или загрузите проект', 'warn'); return; }
  const overlay = $('#progressOverlay');
  overlay.classList.remove('hidden');
  $('#progressText').textContent = 'Расчёт гидродинамики и массопереноса…';
  $('#progressBar').value = 10;
  status('Идёт расчёт…');
  try {
    const body = {
      interpolation: $('#geoForm').interpolation.value,
      variogram_model: $('#geoForm').variogram_model.value,
      variogram_sill: parseFloat($('#geoForm').variogram_sill.value) || 1,
      variogram_range: parseFloat($('#geoForm').variogram_range.value) || 150,
      variogram_nugget: parseFloat($('#geoForm').variogram_nugget.value) || 0,
      anisotropy: [
        parseFloat($('#geoForm').aniso_x.value) || 1,
        parseFloat($('#geoForm').aniso_y.value) || 1,
        parseFloat($('#geoForm').aniso_z.value) || 1,
      ],
      output_every: 30,
      acid_threshold: 5,
      kv_ratio: parseFloat($('#geoForm').kv_ratio.value) || 0.2,
      run_hydro: true, run_flowlines: true, run_chemistry: true, run_reserves: true,
    };
    $('#progressBar').value = 40;
    const res = await api(`/api/projects/${state.pid}/run`, {
      method: 'POST', body: JSON.stringify(body),
    });
    $('#progressBar').value = 80;
    state.results = res;
    await refreshProject();
    await loadFields();
    $('#progressBar').value = 100;
    toast(`Расчёт выполнен за ${res.elapsed_sec ?? '?'} с`);
    status('Расчёт завершён');
  } catch (e) {
    toast(`Ошибка расчёта: ${e.message}`, 'error');
    status('Ошибка расчёта');
    console.error(e);
  } finally {
    overlay.classList.add('hidden');
  }
}

// ------------------------------------------------------------------ поля и срезы

const FIELD_META = {
  head: { title: 'Напор, м', unit: 'м' },
  speed: { title: 'Скорость, м/сут', unit: 'м/сут' },
  k: { title: 'Проницаемость, м/сут', unit: 'м/сут' },
  porosity: { title: 'Пористость', unit: 'д.е.' },
  uranium: { title: 'Уран в руде, кг/т', unit: 'кг/т' },
  acid_final: { title: 'Кислотность раствора, г/л', unit: 'г/л' },
  uranium_final: { title: 'Уран в растворе, мг/л', unit: 'мг/л' },
  residual_final: { title: 'Остаточный уран, кг/т', unit: 'кг/т' },
};

async function loadFields() {
  const name = $('#fieldSelect').value;
  const grid = state.project.grid;
  const layer = parseInt($('#layerSlider').value, 10);
  const key = `${name}:${layer}`;
  let data = state.sliceCache[key];
  if (!data) {
    data = await api(`/api/projects/${state.pid}/field?name=${name}&axis=z&index=${layer}`);
    state.sliceCache[key] = data;
  }
  state.viewer.setSlice(data.data, grid.nx, grid.ny, grid, layer, {
    opacity: parseInt($('#sliceOpacity').value, 10) / 100,
    min: data.min, max: data.max,
  });
  const meta = FIELD_META[name] || { title: name };
  const mn = data.min, mx = data.max;
  $('#legMin').textContent = fmtNum(mn);
  $('#legMid').textContent = fmtNum((mn + mx) / 2);
  $('#legMax').textContent = fmtNum(mx);
  $('#legTitle').textContent = `${meta.title} · слой ${layer + 1}`;
  drawLegend($('#legendCanvas'), mn, mx, (t) => cssColor(t));
  state.sliceStats = { min: mn, max: mx, data: data.data };
  if ($('#rtab-slice').classList.contains('active')) drawSliceMap();
}

function fmtNum(v) {
  const a = Math.abs(v);
  if (a >= 1000) return v.toFixed(0);
  if (a >= 10) return v.toFixed(1);
  if (a >= 1) return v.toFixed(2);
  if (a === 0) return '0';
  return v.toFixed(3);
}

function updateLayerSlider() {
  const grid = state.project.grid;
  const sl = $('#layerSlider');
  sl.max = String(grid.nz - 1);
  // по умолчанию показываем средний слой пласта
  if (!state.layerTouched) {
    sl.value = String(Math.max(0, Math.floor(grid.nz / 2) - 1));
  }
  sl.value = String(Math.min(parseInt(sl.value, 10) || 0, grid.nz - 1));
  $('#layerLabel').textContent = `${parseInt(sl.value, 10) + 1} / ${grid.nz}`;
}

function drawSliceMap() {
  if (!state.sliceStats) return;
  const g = state.project.grid;
  const { min, max, data } = state.sliceStats;
  const layer = parseInt($('#layerSlider').value, 10);
  heatmap($('#sliceCanvas'), data, g.nx, g.ny, {
    norm: (v) => (max > min ? (v - min) / (max - min) : 0.5),
    color: (t) => cssColor(t),
    xMin: `${Math.round(g.x0)} м`,
    xMax: `${Math.round(g.x0 + g.nx * g.dx)} м`,
    yMin: `${Math.round(g.y0)} м`,
    yMax: `${Math.round(g.y0 + g.ny * g.dy)} м`,
    title: `Слой ${layer + 1} (z = ${Math.round(g.z0 + (layer + 0.5) * g.dz)} м)`,
  });
  $('#sliceInfo').textContent =
    `диапазон ${fmtNum(min)} … ${fmtNum(max)}`;
}

// ------------------------------------------------------------------ кривые

function renderCurves() {
  const c = state.results?.curves;
  if (!c) return;
  const t = c.time;
  lineChart($('#curveChart'), [
    { values: c.pr_uranium, color: '#ff5252', x: t, axis: 'left' },
    { values: c.vr_acid, color: '#4d8dff', x: t, axis: 'right' },
  ], { dualAxis: true, xLabel: 'время, сут' });
  lineChart($('#extractChart'), [
    { values: c.cumulative_u, color: '#8bc34a', x: t, axis: 'left' },
  ], { xLabel: 'время, сут' });
}

function renderSummary() {
  const r = state.results;
  const box = $('#summaryBox');
  if (!r || !r.hydro) {
    box.innerHTML = '<p class="muted">Выполните расчёт, чтобы увидеть результаты.</p>';
    return;
  }
  const h = r.hydro;
  const res = r.reserves || {};
  const ps = r.property_stats || {};
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
  const rrows = res.reserve_kg != null ? [
    ['Запас U в закисляемом объёме, т', fmtNum(res.reserve_t)],
    ['Закисляемый объём, м³', fmtNum(res.acidified_volume)],
    ['Площадь, га', fmtNum(res.area_ha)],
    ['Средневзвешенное содержание, кг/т', fmtNum(res.avg_grade)],
    ['Средняя мощность, м', fmtNum(res.avg_thickness)],
    ['Средневзвешенный метропроцент', fmtNum(res.weighted_metre_percent)],
  ] : [];

  box.innerHTML = `
    <h3>Гидродинамика</h3>
    <table class="kv">${rows.map(([k, v]) => `<tr><td>${k}</td><td class="mono">${v}</td></tr>`).join('')}</table>
    ${rrows.length ? `<h3>Запасы</h3><table class="kv">${rrows.map(([k, v]) => `<tr><td>${k}</td><td class="mono">${v}</td></tr>`).join('')}</table>` : ''}
    ${h.balance_ok ? '' : '<p class="warn">⚠ Невязка баланса превышает 1% — проверьте расходы скважин и граничные условия.</p>'}`;
}

// ------------------------------------------------------------------ линии тока

async function loadStreamlines() {
  try {
    const res = await api(`/api/projects/${state.pid}/streamlines?limit=400`);
    state.flowlines = res.streamlines || [];
    state.viewer.setStreamlines(state.flowlines, state.project.grid);
    const sel = $('#flowlineSelect');
    sel.innerHTML = '';
    state.flowlines.slice(0, 200).forEach((l, i) => {
      const o = document.createElement('option');
      o.value = String(i);
      o.textContent = `${l.well_in || '—'} → ${l.well_out || 'выход'} · ${Math.round(l.residence)} сут`;
      sel.appendChild(o);
    });
  } catch (e) {
    console.warn('Нет линий тока', e);
  }
}

async function selectFlowline() {
  const idx = parseInt($('#flowlineSelect').value || '0', 10);
  if (!state.project) return;
  try {
    const res = await api(`/api/projects/${state.pid}/flowline-chemistry?index=${idx}`);
    lineChart($('#flowlineChart'), [
      { values: res.uranium, color: '#ff5252', x: res.times, axis: 'left' },
      { values: res.acid, color: '#4d8dff', x: res.times, axis: 'left' },
      { values: res.rel.map((v) => v * 100), color: '#8bc34a', x: res.times, axis: 'left' },
    ], { xLabel: 'время пребывания, сут' });
    const wname = (id) => {
      const w = state.project.wells.find((x) => x.id === id);
      return w ? w.name : (id || '—');
    };
    $('#flowlineInfo').innerHTML = `
      <div>Линия: ${wname(res.well_in)} → ${wname(res.well_out) === '—' ? 'выход из блока' : wname(res.well_out)}</div>
      <div>Время пребывания: ${fmtNum(res.residence)} сут</div>
      <div>Кислота на выходе: ${fmtNum(res.acid[res.acid.length - 1])} г/л</div>
      <div>Уран в растворе: ${fmtNum(res.uranium[res.uranium.length - 1])} мг/л</div>
      <div>Остаток U в руде: ${(100 * (res.rel[res.rel.length - 1])).toFixed(1)} %</div>`;
  } catch (e) {
    $('#flowlineInfo').textContent = 'Не удалось рассчитать линию тока';
  }
}

// ------------------------------------------------------------------ импорт

function initImport() {
  $('#btnImport').addEventListener('click', () => $('#importModal').classList.remove('hidden'));
  $('#btnCloseImport').addEventListener('click', () => $('#importModal').classList.add('hidden'));
  $('#importFile').addEventListener('change', async (e) => {
    const file = e.target.files[0];
    if (!file) return;
    const fd = new FormData();
    fd.append('file', file);
    status('Загрузка файла…');
    try {
      const res = await fetch(`/api/projects/${state.pid}/import`, { method: 'POST', body: fd });
      if (!res.ok) throw new Error((await res.json()).detail || 'ошибка загрузки');
      const data = await res.json();
      state.importData = data;
      renderMapping(data);
      $('#importCount').textContent = `Распознано скважин: ${data.count}`;
      status('Файл загружен');
    } catch (err) {
      toast(`Ошибка импорта: ${err.message}`, 'error');
    }
  });
  $('#btnApplyImport').addEventListener('click', async () => {
    if (!state.importData) return;
    try {
      await api(`/api/projects/${state.pid}/import/apply`, {
        method: 'POST',
        body: JSON.stringify({
          rows: state.importData.rows,
          mapping: state.importData.mapping,
          as_properties: true,
        }),
      });
      await refreshProject();
      $('#importModal').classList.add('hidden');
      toast('Данные загружены в проект');
    } catch (err) {
      toast(`Ошибка: ${err.message}`, 'error');
    }
  });
}

function renderMapping(data) {
  const box = $('#mappingRows');
  box.innerHTML = '';
  const labels = {
    name: 'Наименование', x: 'X', y: 'Y', z: 'Z (отметка)', k: 'Проницаемость',
    porosity: 'Пористость', uranium: 'Уран', clay: 'Глина', lithology: 'Литология',
    rate: 'Дебит', well_type: 'Тип', acid_in: 'Кислотность ВР',
    screen_top: 'Верх фильтра', screen_bottom: 'Низ фильтра',
  };
  const cols = ['', ...data.columns];
  for (const [field, label] of Object.entries(labels)) {
    const row = document.createElement('div');
    row.className = 'map-row';
    row.innerHTML = `<span>${label}</span><select data-field="${field}">${
      cols.map((c) => `<option value="${c}"${c === data.mapping[field] ? ' selected' : ''}>${c || '—'}</option>`).join('')
    }</select>`;
    box.appendChild(row);
  }
  $('#importMapping').classList.remove('hidden');
  box.addEventListener('change', () => {
    const mapping = {};
    box.querySelectorAll('select').forEach((s) => { mapping[s.dataset.field] = s.value || null; });
    state.importData.mapping = mapping;
    // пересчёт строк с новым сопоставлением
    state.importData.rows = remapRows(state.importData, mapping);
  });
}

function remapRows(data, mapping) {
  const out = [];
  for (const rec of data.rows) {
    const r = { ...rec };
    for (const [field, col] of Object.entries(mapping)) {
      if (!col) continue;
      const idx = data.columns.indexOf(col);
      if (idx >= 0 && data.preview) { /* preview неполон — сохраняем исходное значение */ }
    }
    out.push(r);
  }
  return out;
}

// ------------------------------------------------------------------ сохранение

async function saveProject() {
  if (!state.pid) return;
  try {
    const r = await api(`/api/projects/${state.pid}/save`, { method: 'POST' });
    toast(`Сохранено: ${r.path}`);
  } catch (e) {
    toast(`Ошибка сохранения: ${e.message}`, 'error');
  }
}

async function openProject() {
  const input = document.createElement('input');
  input.type = 'file';
  input.accept = '.s3dp';
  input.addEventListener('change', async () => {
    const file = input.files[0];
    if (!file) return;
    const fd = new FormData();
    fd.append('file', file);
    try {
      const p = await fetch('/api/projects/open', { method: 'POST', body: fd });
      if (!p.ok) throw new Error((await p.json()).detail || 'ошибка');
      const data = await p.json();
      state.pid = data.id;
      await refreshProject();
      toast('Проект загружен');
    } catch (e) {
      toast(`Ошибка: ${e.message}`, 'error');
    }
  });
  input.click();
}

// ------------------------------------------------------------------ формы

function initForms() {
  $('#wellForm').addEventListener('submit', async (e) => {
    e.preventDefault();
    const f = e.target;
    const wells = state.project.wells.map((w) => ({
      id: w.id, name: w.name, well_type: w.well_type, x: w.x, y: w.y, z: w.z ?? 0,
      screen_top: w.screen_top, screen_bottom: w.screen_bottom,
      rate: w.rate, acid_in: w.acid_in, screens: w.screens || [],
    }));
    const idx = state.editingWell
      ? wells.findIndex((w) => w.id === state.editingWell.id) : -1;
    const item = {
      id: state.editingWell ? state.editingWell.id : undefined,
      name: f.name.value, well_type: f.well_type.value,
      x: parseFloat(f.x.value), y: parseFloat(f.y.value), z: parseFloat(f.z.value),
      acid_in: parseFloat(f.acid_in.value),
      screen_top: f.screen_top.value === '' ? null : parseFloat(f.screen_top.value),
      screen_bottom: f.screen_bottom.value === '' ? null : parseFloat(f.screen_bottom.value),
      rate: parseFloat(f.rate.value),
    };
    if (idx >= 0) wells[idx] = { ...wells[idx], ...item }; else wells.push(item);
    await api(`/api/projects/${state.pid}/wells`, {
      method: 'PUT', body: JSON.stringify(wells),
    });
    state.editingWell = null;
    f.classList.add('hidden');
    await refreshProject();
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

  $('#gridForm').addEventListener('submit', async (e) => {
    e.preventDefault();
    const f = e.target;
    const g = {
      x0: parseFloat(f.x0.value), y0: parseFloat(f.y0.value), z0: parseFloat(f.z0.value),
      dx: parseFloat(f.dx.value), dy: parseFloat(f.dy.value), dz: parseFloat(f.dz.value),
      nx: parseInt(f.nx.value), ny: parseInt(f.ny.value), nz: parseInt(f.nz.value),
    };
    await api(`/api/projects/${state.pid}/grid`, { method: 'PUT', body: JSON.stringify(g) });
    await refreshProject();
    toast('Сетка обновлена');
  });

  $('#chemForm').addEventListener('submit', async (e) => {
    e.preventDefault();
    const f = e.target;
    const body = { ...state.project.chem };
    ['acid_concentration', 'fe3_conc', 'feoh3_conc', 'caco3_conc', 'integration_days',
      'temperature', 'activation_energy', 'beta', 'k1', 'k2', 'k3', 'k4',
      'flowline_k1', 'flowline_k2', 'flowline_k3', 'flowline_k4', 'flowline_beta',
      'flowline_acid'].forEach((k) => {
      if (f[k] && f[k].value !== '') body[k] = parseFloat(f[k].value);
    });
    await api(`/api/projects/${state.pid}/chemistry`, { method: 'PUT', body: JSON.stringify(body) });
    await refreshProject();
    toast('Параметры кинетики сохранены');
  });

  $('#bcForm').addEventListener('submit', async (e) => {
    e.preventDefault();
    const faces = {}, heads = {};
    $$('#bcRows select').forEach((s) => { faces[s.dataset.face] = s.value; });
    $$('#bcRows input').forEach((i) => { heads[i.dataset.head] = parseFloat(i.value) || 0; });
    await api(`/api/projects/${state.pid}/bc`, {
      method: 'PUT',
      body: JSON.stringify({ boundary: 'neumann', boundary_head: 0, bc_faces: faces, bc_heads: heads }),
    });
    await refreshProject();
    toast('Граничные условия применены');
  });

  $('#geoForm').addEventListener('submit', (e) => {
    e.preventDefault();
    toast('Настройки интерполяции сохранены и будут применены при расчёте');
  });
}

// ------------------------------------------------------------------ init

function initViewerControls() {
  $('#fieldSelect').addEventListener('change', () => { if (state.pid) loadFields(); });
  $('#layerSlider').addEventListener('input', () => {
    state.layerTouched = true;
    $('#layerLabel').textContent = `${parseInt($('#layerSlider').value, 10) + 1} / ${state.project.grid.nz}`;
    if (state.pid) loadFields();
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

async function init() {
  initTabs();
  state.viewer = new Viewer3D($('#viewer'));
  initForms();
  initImport();
  initViewerControls();

  $('#btnDemo').addEventListener('click', loadDemo);
  $('#btnNew').addEventListener('click', newProject);
  $('#btnRun').addEventListener('click', runCalc);
  $('#btnSave').addEventListener('click', saveProject);
  $('#btnOpen').addEventListener('click', openProject);

  // Первичная загрузка демо-проекта
  await loadDemo();
  const res = await api(`/api/projects/${state.pid}/run`, {
    method: 'POST',
    body: JSON.stringify({ run_hydro: true, run_flowlines: true, run_chemistry: true, run_reserves: true }),
  });
  state.results = res;
  await refreshProject();
  await loadFields();
  await loadStreamlines();
  renderCurves();
  selectFlowline();
  status('Демо-блок рассчитан');
}

init().catch((e) => {
  console.error(e);
  status(`Ошибка инициализации: ${e.message}`);
  toast(`Ошибка: ${e.message}`, 'error');
});
