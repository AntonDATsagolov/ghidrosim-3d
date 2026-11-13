/**
 * Вычислительное ядро «ГидроСимулятор 3D» — порт на JavaScript.
 *
 * Полностью совпадает с Python-реализацией (backend/app): те же уравнения,
 * те же конечные объёмы, та же схема адвекции-дисперсии и те же кинетические
 * уравнения. Модуль не зависит от DOM — его можно запускать и в браузере,
 * и в Node.js (для сверки с эталоном).
 *
 * Плоские массивы хранят поля в порядке C: index = i + nx*(j + ny*k),
 * где i — индекс по X, j — по Y, k — по Z.
 */

export const ORE_DENSITY = 2200;     // кг/м3
export const ACID_PER_U = 2.0;       // кг H2SO4 на кг U
export const FE3_PER_U = 0.47;       // кг Fe(III) на кг U
export const T0_REF = 293.15;        // К
export const R_GAS = 8.314;          // Дж/(моль·К)

// ------------------------------------------------------------------ утилиты

/**
 * Плоский индекс ячейки в порядке C (как numpy): offset = i*(ny*nz) + j*nz + k.
 * Именно такой порядок используется во всех полях и на гранях, поэтому
 * значения совпадают с эталонной Python-реализацией.
 */
export function idx3(i, j, k, ny, nz) { return i * ny * nz + j * nz + k; }

export function clamp(v, lo, hi) { return v < lo ? lo : (v > hi ? hi : v); }

export function mean(arr) {
  let s = 0;
  for (let i = 0; i < arr.length; i++) s += arr[i];
  return arr.length ? s / arr.length : 0;
}

export function percentile(arr, p) {
  if (!arr.length) return 0;
  const a = Array.from(arr).sort((x, y) => x - y);
  const idx = clamp(p / 100 * (a.length - 1), 0, a.length - 1);
  const lo = Math.floor(idx), hi = Math.ceil(idx);
  return a[lo] + (a[hi] - a[lo]) * (idx - lo);
}

/** Решение плотной симметричной системы методом Холецкого (n ≤ 40). */
function cholSolve(A, b, n) {
  const L = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let s = A[i * n + j];
      for (let k = 0; k < j; k++) s -= L[i * n + k] * L[j * n + k];
      if (i === j) {
        if (s <= 1e-14) s = 1e-14;
        L[i * n + j] = Math.sqrt(s);
      } else {
        L[i * n + j] = s / L[j * n + j];
      }
    }
  }
  const y = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let s = b[i];
    for (let k = 0; k < i; k++) s -= L[i * n + k] * y[k];
    y[i] = s / L[i * n + i];
  }
  const x = new Float64Array(n);
  for (let i = n - 1; i >= 0; i--) {
    let s = y[i];
    for (let k = i + 1; k < n; k++) s -= L[k * n + i] * x[k];
    x[i] = s / L[i * n + i];
  }
  return x;
}

/**
 * Метод Гаусса с частичным выбором ведущего элемента.
 * Используется для расширенной системы кригинга, которая не является
 * положительно определённой (в углу 0), поэтому Холецкий неприменим.
 */
function luSolve(A, b, n) {
  const M = Float64Array.from(A);
  const y = Float64Array.from(b);
  for (let col = 0; col < n; col++) {
    // ведущий элемент
    let piv = col, best = Math.abs(M[col * n + col]);
    for (let r = col + 1; r < n; r++) {
      const v = Math.abs(M[r * n + col]);
      if (v > best) { best = v; piv = r; }
    }
    if (best < 1e-300) continue;
    if (piv !== col) {
      for (let c = 0; c < n; c++) {
        const t = M[col * n + c]; M[col * n + c] = M[piv * n + c]; M[piv * n + c] = t;
      }
      const t = y[col]; y[col] = y[piv]; y[piv] = t;
    }
    const d = M[col * n + col];
    for (let r = col + 1; r < n; r++) {
      const f = M[r * n + col] / d;
      if (f === 0) continue;
      for (let c = col; c < n; c++) M[r * n + c] -= f * M[col * n + c];
      y[r] -= f * y[col];
    }
  }
  const x = new Float64Array(n);
  for (let i = n - 1; i >= 0; i--) {
    let s = y[i];
    for (let k = i + 1; k < n; k++) s -= M[i * n + k] * x[k];
    const d = M[i * n + i];
    x[i] = Math.abs(d) < 1e-300 ? 0 : s / d;
  }
  return x;
}

/** Метод наименьших квадратов для недоопределённой системы (нормальные уравнения). */
function lstsqSolve(A, b, m, n) {
  const AtA = new Float64Array(n * n);
  const Atb = new Float64Array(n);
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      const a = A[i * n + j];
      Atb[j] += a * b[i];
      for (let k = j; k < n; k++) AtA[j * n + k] += a * A[i * n + k];
    }
  }
  for (let j = 0; j < n; j++) for (let k = 0; k < j; k++) AtA[j * n + k] = AtA[k * n + j];
  return cholSolve(AtA, Atb, n);
}

// ------------------------------------------------------------------ вариограмма

export function variogram(model, sill, range, nugget, h) {
  const r = Math.max(range, 1e-9);
  let f;
  switch (model) {
    case 'exponential': f = 1 - Math.exp(-3 * h / r); break;
    case 'gaussian': f = 1 - Math.exp(-3 * (h / r) ** 2); break;
    case 'linear': f = clamp(h / r, 0, 1); break;
    default: { const hr = clamp(h / r, 0, 1); f = 1.5 * hr - 0.5 * hr ** 3; }
  }
  return nugget + sill * f;
}

// ------------------------------------------------------------------ интерполяция

function applyAnisotropy(dx, dy, dz, aniso) {
  return Math.sqrt((dx / aniso[0]) ** 2 + (dy / aniso[1]) ** 2 + (dz / aniso[2]) ** 2);
}

/** Обратные взвешенные расстояния. */
export function idw(points, values, targets, opts) {
  const power = opts.idwPower ?? 2;
  const aniso = opts.anisotropy ?? [1, 1, 1];
  const out = new Float64Array(targets.length);
  for (let t = 0; t < targets.length; t++) {
    const d = new Float64Array(points.length);
    for (let p = 0; p < points.length; p++) {
      d[p] = applyAnisotropy(
        targets[t][0] - points[p][0], targets[t][1] - points[p][1],
        targets[t][2] - points[p][2], aniso);
    }
    const order = Array.from(d.keys()).sort((a, b) => d[a] - d[b]);
    const k = Math.min(opts.maxNeighbors ?? 12, order.length);
    let sw = 0, sv = 0;
    for (let i = 0; i < k; i++) {
      const j = order[i];
      const w = 1 / Math.max(d[j], 1e-6) ** power;
      sw += w; sv += w * values[j];
    }
    out[t] = sw > 0 ? sv / sw : 0;
  }
  return out;
}

/** Ординарный кригинг (система в форме вариограммы). */
export function kriging(points, values, targets, opts) {
  const vg = opts.variogram ?? { model: 'spherical', sill: 1, range: 150, nugget: 0 };
  const aniso = opts.anisotropy ?? [1, 1, 1];
  const n = points.length;
  const maxPts = 16;
  const out = new Float64Array(targets.length);

  // Матрица вариограмм между исходными точками
  const K = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const h = applyAnisotropy(
        points[i][0] - points[j][0], points[i][1] - points[j][1],
        points[i][2] - points[j][2], aniso);
      K[i * n + j] = variogram(vg.model, vg.sill, vg.range, vg.nugget, h);
    }
  }

  // Ближайшие соседи для каждой цели
  const neigh = [];
  for (let t = 0; t < targets.length; t++) {
    const d = new Float64Array(n);
    for (let p = 0; p < n; p++) {
      d[p] = applyAnisotropy(
        targets[t][0] - points[p][0], targets[t][1] - points[p][1],
        targets[t][2] - points[p][2], aniso);
    }
    const order = Array.from(d.keys()).sort((a, b) => d[a] - d[b]);
    neigh.push(order.slice(0, Math.min(maxPts, n)));
  }

  for (let t = 0; t < targets.length; t++) {
    const idx = neigh[t];
    const m = idx.length;
    const A = new Float64Array((m + 1) * (m + 1));
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < m; j++) A[i * (m + 1) + j] = K[idx[i] * n + idx[j]];
      A[i * (m + 1) + m] = 1;
      A[m * (m + 1) + i] = 1;
    }
    // регуляризация диагонали
    let tr = 0;
    for (let i = 0; i < m; i++) tr += A[i * (m + 1) + i];
    const nug = Math.max(vg.nugget, 1e-8) + Math.max(Math.abs(tr) / m, 1e-8) * 1e-6;
    for (let i = 0; i < m; i++) A[i * (m + 1) + i] += nug;

    const b = new Float64Array(m + 1);
    for (let i = 0; i < m; i++) {
      const h = applyAnisotropy(
        targets[t][0] - points[idx[i]][0], targets[t][1] - points[idx[i]][1],
        targets[t][2] - points[idx[i]][2], aniso);
      b[i] = variogram(vg.model, vg.sill, vg.range, vg.nugget, h);
    }
    b[m] = 1;

    let w;
    try {
      w = luSolve(A, b, m + 1);
    } catch (e) {
      w = lstsqSolve(A, b, m + 1, m + 1);
    }
    if (!Number.isFinite(w[0])) w = lstsqSolve(A, b, m + 1, m + 1);
    let est = 0;
    for (let i = 0; i < m; i++) est += w[i] * values[idx[i]];
    out[t] = est;
  }
  return out;
}

export function interpolate(points, values, targets, opts) {
  if (!points.length) return new Float64Array(targets.length);
  const vals = [];
  const pts = [];
  for (let i = 0; i < points.length; i++) {
    if (Number.isFinite(values[i])) { pts.push(points[i]); vals.push(values[i]); }
  }
  if (!pts.length) return new Float64Array(targets.length);
  return opts.method === 'idw'
    ? idw(pts, vals, targets, opts)
    : kriging(pts, vals, targets, opts);
}

// ------------------------------------------------------------------ сетка

export function makeGrid(spec) {
  const g = {
    x0: spec.x0 ?? 0, y0: spec.y0 ?? 0, z0: spec.z0 ?? 0,
    nx: spec.nx, ny: spec.ny, nz: spec.nz,
    dx: spec.dx, dy: spec.dy, dz: spec.dz,
  };
  g.n = g.nx * g.ny * g.nz;
  g.cellVolume = g.dx * g.dy * g.dz;
  g.bounds = {
    x_min: g.x0, x_max: g.x0 + g.nx * g.dx,
    y_min: g.y0, y_max: g.y0 + g.ny * g.dy,
    z_min: g.z0, z_max: g.z0 + g.nz * g.dz,
  };
  return g;
}

export function gridCenters(g) {
  const x = new Float64Array(g.nx), y = new Float64Array(g.ny), z = new Float64Array(g.nz);
  for (let i = 0; i < g.nx; i++) x[i] = g.x0 + (i + 0.5) * g.dx;
  for (let j = 0; j < g.ny; j++) y[j] = g.y0 + (j + 0.5) * g.dy;
  for (let k = 0; k < g.nz; k++) z[k] = g.z0 + (k + 0.5) * g.dz;
  return { x, y, z };
}

/** Целевые точки для интерполяции (по одной на ячейку). */
export function gridTargets(g) {
  const c = gridCenters(g);
  const out = new Array(g.nx * g.ny * g.nz);
  for (let k = 0; k < g.nz; k++)
    for (let j = 0; j < g.ny; j++)
      for (let i = 0; i < g.nx; i++)
        out[idx3(i, j, k, g.ny, g.nz)] = [c.x[i], c.y[j], c.z[k]];
  return out;
}

// ------------------------------------------------------------------ свойства

export function buildPropertyModels(project, opts) {
  const g = project.grid;
  const nCells = g.n ?? g.nx * g.ny * g.nz;
  const targets = gridTargets(g);
  // Опции принимаются и в camelCase, и в snake_case
  const pick = (...names) => {
    for (const nm of names) if (opts[nm] !== undefined) return opts[nm];
    return undefined;
  };
  const settings = {
    method: pick('interpolation', 'interpolationMethod') ?? 'kriging',
    variogram: {
      model: pick('variogramModel', 'variogram_model') ?? 'spherical',
      sill: pick('variogramSill', 'variogram_sill') ?? 1,
      range: pick('variogramRange', 'variogram_range') ?? 150,
      nugget: pick('variogramNugget', 'variogram_nugget') ?? 0,
    },
    anisotropy: pick('anisotropy') ?? [1, 1, 1],
    maxNeighbors: 12,
  };

  const interp = (list, def) => {
    if (!list.length) return new Float64Array(nCells).fill(def);
    const pts = list.map((r) => [r[0], r[1], r[2]]);
    const vals = list.map((r) => r[3]);
    const res = interpolate(pts, vals, targets, settings);
    const lo = percentile(vals, 1), hi = percentile(vals, 99);
    const a = Math.min(lo, hi), b = Math.max(lo, hi);
    for (let i = 0; i < nCells; i++) {
      let v = Number.isFinite(res[i]) ? res[i] : def;
      if (b > a) v = clamp(v, a, b);
      res[i] = v;
    }
    return res;
  };

  const lithology = new Array(nCells).fill('');
  const lithData = project.props.lithology || [];
  if (lithData.length) {
    const pts = lithData.map((r) => [r[0], r[1], r[2]]);
    const labels = lithData.map((r) => r[3]);
    for (let i = 0; i < nCells; i++) {
      const t = targets[i];
      let best = 0, bd = Infinity;
      for (let p = 0; p < pts.length; p++) {
        const d = (pts[p][0] - t[0]) ** 2 + (pts[p][1] - t[1]) ** 2 + (pts[p][2] - t[2]) ** 2;
        if (d < bd) { bd = d; best = p; }
      }
      lithology[i] = labels[best];
    }
  }

  return {
    k: interp(project.props.permeability, 1.0),
    porosity: interp(project.props.porosity, 0.15),
    uranium: interp(project.props.uranium, 0.0),
    clay: interp(project.props.clay, 0.1),
    lithology,
  };
}

// ------------------------------------------------------------------ гидродинамика

const FACES = {
  x_min: [0, 0], x_max: [0, -1],
  y_min: [1, 0], y_max: [1, -1],
  z_min: [2, 0], z_max: [2, -1],
};

export function wellScreenCells(well, g) {
  const out = [];
  const c = gridCenters(g);
  const intervals = well.screens && well.screens.length
    ? well.screens
    : (well.screen ? [well.screen] : []);
  if (!well.trajectory || !well.trajectory.length || !intervals.length) return out;
  const xh = well.trajectory[0].x, yh = well.trajectory[0].y;
  for (const iv of intervals) {
    for (let k = 0; k < g.nz; k++) {
      if (c.z[k] < iv.top - 1e-9 || c.z[k] > iv.bottom + 1e-9) continue;
      for (let j = 0; j < g.ny; j++) {
        if (Math.abs(c.y[j] - yh) > 0.5 * g.dy + 1e-9) continue;
        for (let i = 0; i < g.nx; i++) {
          if (Math.abs(c.x[i] - xh) > 0.5 * g.dx + 1e-9) continue;
          out.push(idx3(i, j, k, g.ny, g.nz));
        }
      }
    }
  }
  return out;
}

/** Источники по ячейкам (м3/сут) и карта «скважина → ячейки». */
export function buildSources(g, wells, injFactor = 1, prodFactor = 1) {
  const q = new Float64Array(g.n);
  const cellsMap = {};
  for (const w of wells) {
    const cells = wellScreenCells(w, g);
    cellsMap[w.id] = cells;
    if (!cells.length || Math.abs(w.rate) < 1e-12) continue;
    let total;
    if (w.well_type === 'injection') total = w.rate * injFactor;
    else if (w.well_type === 'production') total = -Math.abs(w.rate) * prodFactor;
    else continue;
    const per = total / cells.length;
    for (const c of cells) q[c] += per;
  }
  return { q, cellsMap };
}

/** Проводимости на гранях: T = k_гарм * A / Δx (м3/сут на единицу напора). */
function faceConductance(g, k, axis, kvRatio) {
  const out = new Float64Array(
    axis === 0 ? (g.nx - 1) * g.ny * g.nz
      : axis === 1 ? g.nx * (g.ny - 1) * g.nz
        : g.nx * g.ny * (g.nz - 1));
  const harm = (a, b) => {
    const s = a + b;
    return s > 0 ? 2 * a * b / s : 0;
  };
  if (axis === 0) {
    const area = g.dy * g.dz, dist = g.dx;
    for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny; jj++)
      for (let ii = 0; ii < g.nx - 1; ii++)
        out[ii * (g.ny * g.nz) + jj * g.nz + kk] = harm(
          k[idx3(ii, jj, kk, g.ny, g.nz)], k[idx3(ii + 1, jj, kk, g.ny, g.nz)]) * area / dist;
  } else if (axis === 1) {
    const area = g.dx * g.dz, dist = g.dy;
    for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny - 1; jj++)
      for (let ii = 0; ii < g.nx; ii++)
        out[ii * ((g.ny - 1) * g.nz) + jj * g.nz + kk] = harm(
          k[idx3(ii, jj, kk, g.ny, g.nz)], k[idx3(ii, jj + 1, kk, g.ny, g.nz)]) * area / dist;
  } else {
    const area = g.dx * g.dy, dist = g.dz;
    for (let kk = 0; kk < g.nz - 1; kk++) for (let jj = 0; jj < g.ny; jj++)
      for (let ii = 0; ii < g.nx; ii++)
        out[ii * (g.ny * (g.nz - 1)) + jj * (g.nz - 1) + kk] = harm(
          k[idx3(ii, jj, kk, g.ny, g.nz)], k[idx3(ii, jj, kk + 1, g.ny, g.nz)]) * kvRatio * area / dist;
  }
  return out;
}

/** Проницаемость на грани между соседними ячейками (гармоническое среднее). */
function facePerm(g, k, axis) {
  const harm = (a, b) => {
    const s = a + b;
    return s > 0 ? 2 * a * b / s : 0;
  };
  const size = axis === 0 ? (g.nx - 1) * g.ny * g.nz
    : axis === 1 ? g.nx * (g.ny - 1) * g.nz
      : g.nx * g.ny * (g.nz - 1);
  const out = new Float64Array(size);
  if (axis === 0) {
    for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny; jj++)
      for (let ii = 0; ii < g.nx - 1; ii++)
        out[ii * (g.ny * g.nz) + jj * g.nz + kk] = harm(
          k[idx3(ii, jj, kk, g.ny, g.nz)], k[idx3(ii + 1, jj, kk, g.ny, g.nz)]);
  } else if (axis === 1) {
    for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny - 1; jj++)
      for (let ii = 0; ii < g.nx; ii++)
        out[ii * ((g.ny - 1) * g.nz) + jj * g.nz + kk] = harm(
          k[idx3(ii, jj, kk, g.ny, g.nz)], k[idx3(ii, jj + 1, kk, g.ny, g.nz)]);
  } else {
    for (let kk = 0; kk < g.nz - 1; kk++) for (let jj = 0; jj < g.ny; jj++)
      for (let ii = 0; ii < g.nx; ii++)
        out[ii * (g.ny * (g.nz - 1)) + jj * (g.nz - 1) + kk] = harm(
          k[idx3(ii, jj, kk, g.ny, g.nz)], k[idx3(ii, jj, kk + 1, g.ny, g.nz)]);
  }
  return out;
}

/**
 * Сборка матрицы оператора (для сверки с эталоном и отладки).
 * Возвращает диагональ и проводимости на гранях в порядке C.
 */
export function buildOperatorMatrix(g, kArr, kvRatio = 1.0) {
  const n = g.nx * g.ny * g.nz;
  const k = new Float64Array(n);
  for (let i = 0; i < n; i++) k[i] = (Number.isFinite(kArr[i]) && kArr[i] > 0) ? kArr[i] : 1e-6;
  const Tx = faceConductance(g, k, 0, kvRatio);
  const Ty = faceConductance(g, k, 1, kvRatio);
  const Tz = faceConductance(g, k, 2, kvRatio);
  const diag = new Float64Array(n);
  for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx; ii++) {
    const c = idx3(ii, jj, kk, g.ny, g.nz);
    if (ii > 0) diag[c] += Tx[(ii - 1) * (g.ny * g.nz) + jj * g.nz + kk];
    if (ii < g.nx - 1) diag[c] += Tx[ii * (g.ny * g.nz) + jj * g.nz + kk];
    if (jj > 0) diag[c] += Ty[ii * ((g.ny - 1) * g.nz) + (jj - 1) * g.nz + kk];
    if (jj < g.ny - 1) diag[c] += Ty[ii * ((g.ny - 1) * g.nz) + jj * g.nz + kk];
    if (kk > 0) diag[c] += Tz[ii * (g.ny * (g.nz - 1)) + jj * (g.nz - 1) + (kk - 1)];
    if (kk < g.nz - 1) diag[c] += Tz[ii * (g.ny * (g.nz - 1)) + jj * (g.nz - 1) + kk];
  }
  return { diag, Tx, Ty, Tz };
}

/**
 * Стационарная фильтрация Дарси в 3D.
 * Уравнение: div(k grad h) + q = 0.
 * Решается методом сопряжённых градиентов с предобусловливанием Жаку.
 */
export function solveDarcy(g, kArr, sources, porosity, bcFaces = {}, bcHeads = {}, kvRatio = 1.0) {
  const n = g.n ?? g.nx * g.ny * g.nz;
  const k = new Float64Array(n);
  for (let i = 0; i < n; i++) k[i] = (Number.isFinite(kArr[i]) && kArr[i] > 0) ? kArr[i] : 1e-6;

  const Tx = faceConductance(g, k, 0, kvRatio);
  const Ty = faceConductance(g, k, 1, kvRatio);
  const Tz = faceConductance(g, k, 2, kvRatio);

  const diag = new Float64Array(n);
  const b = new Float64Array(n);
  for (let i = 0; i < n; i++) b[i] = sources[i];

  const addDiag = (cell, t) => { diag[cell] += t; };

  // Диагональ графа Лапласа
  for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx; ii++) {
    const c = idx3(ii, jj, kk, g.ny, g.nz);
    if (ii > 0) addDiag(c, Tx[(ii - 1) * (g.ny * g.nz) + jj * g.nz + kk]);
    if (ii < g.nx - 1) addDiag(c, Tx[ii * (g.ny * g.nz) + jj * g.nz + kk]);
    if (jj > 0) addDiag(c, Ty[ii * ((g.ny - 1) * g.nz) + (jj - 1) * g.nz + kk]);
    if (jj < g.ny - 1) addDiag(c, Ty[ii * ((g.ny - 1) * g.nz) + jj * g.nz + kk]);
    if (kk > 0) addDiag(c, Tz[ii * (g.ny * (g.nz - 1)) + jj * (g.nz - 1) + (kk - 1)]);
    if (kk < g.nz - 1) addDiag(c, Tz[ii * (g.ny * (g.nz - 1)) + jj * (g.nz - 1) + kk]);
  }

  // Граничные условия Дирихле
  let hasDirichlet = false;
  for (const [face, cond] of Object.entries(bcFaces)) {
    if (cond !== 'dirichlet') continue;
    hasDirichlet = true;
    const [axis, side] = FACES[face];
    const hd = bcHeads[face] ?? 0;
    const half = axis === 0 ? 0.5 * g.dx : axis === 1 ? 0.5 * g.dy : 0.5 * g.dz;
    const area = axis === 0 ? g.dy * g.dz : axis === 1 ? g.dx * g.dz : g.dx * g.dy;
    for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx; ii++) {
      let onFace = false;
      if (axis === 0) onFace = (side === 0 ? ii === 0 : ii === g.nx - 1);
      else if (axis === 1) onFace = (side === 0 ? jj === 0 : jj === g.ny - 1);
      else onFace = (side === 0 ? kk === 0 : kk === g.nz - 1);
      if (!onFace) continue;
      const c = idx3(ii, jj, kk, g.ny, g.nz);
      const kEff = k[c] * (axis === 2 ? kvRatio : 1);
      const t = kEff * area / half;
      diag[c] += t;
      b[c] += t * hd;
    }
  }

  // Чистый Нейман: система согласована только при нулевом балансе.
  // Избыток расхода равномерно уходит через внешние грани.
  const rankOne = new Float64Array(n);   // поправка (alpha/n)·11ᵀ
  const diagPC = new Float64Array(n);    // диагональ предобусловливателя
  let alpha = 0;
  if (!hasDirichlet) {
    let net = 0;
    for (let i = 0; i < n; i++) net += b[i];
    if (Math.abs(net) > 1e-9) {
      const weight = new Float64Array(n);
      for (const [face, [axis, side]] of Object.entries(FACES)) {
        for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx; ii++) {
          let onFace = false;
          if (axis === 0) onFace = (side === 0 ? ii === 0 : ii === g.nx - 1);
          else if (axis === 1) onFace = (side === 0 ? jj === 0 : jj === g.ny - 1);
          else onFace = (side === 0 ? kk === 0 : kk === g.nz - 1);
          if (onFace) weight[idx3(ii, jj, kk, g.ny, g.nz)] += 1;
        }
      }
      let totW = 0;
      for (let i = 0; i < n; i++) totW += weight[i];
      if (totW > 0) for (let i = 0; i < n; i++) b[i] -= net * (weight[i] / totW);
    }
    // Ранговая поправка (alpha/n)·11ᵀ устраняет вырожденность по константе.
    // Решение при этом не искажается: правую часть уже сбалансировали,
    // поэтому искомое решение имеет нулевое среднее и удовлетворяет
    // как A·x = b, так и (A + alpha/n·11ᵀ)·x = b.
    let maxD = 0;
    for (let i = 0; i < n; i++) maxD = Math.max(maxD, Math.abs(diag[i]));
    alpha = Math.max(maxD, 1) * 0.1;
    const c = alpha / n;
    for (let i = 0; i < n; i++) rankOne[i] = c;
    // диагональ предобусловливателя (оператор остаётся без поправки —
    // иначе ранговый член учитывался бы дважды)
    diagPC.set(diag);
    for (let i = 0; i < n; i++) diagPC[i] += c;
  }

  // Оператор A: (A h)_i = diag_i h_i - sum_j T_ij h_j (+ ранговая поправка)
  const applyOp = (h, out) => {
    let sum = 0;
    if (alpha > 0) { for (let i = 0; i < n; i++) sum += h[i]; sum *= alpha / n; }
    for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx; ii++) {
      const c = idx3(ii, jj, kk, g.ny, g.nz);
      let v = diag[c] * h[c];
      if (ii > 0) v -= Tx[(ii - 1) * (g.ny * g.nz) + jj * g.nz + kk] * h[idx3(ii - 1, jj, kk, g.ny, g.nz)];
      if (ii < g.nx - 1) v -= Tx[ii * (g.ny * g.nz) + jj * g.nz + kk] * h[idx3(ii + 1, jj, kk, g.ny, g.nz)];
      if (jj > 0) v -= Ty[ii * ((g.ny - 1) * g.nz) + (jj - 1) * g.nz + kk] * h[idx3(ii, jj - 1, kk, g.ny, g.nz)];
      if (jj < g.ny - 1) v -= Ty[ii * ((g.ny - 1) * g.nz) + jj * g.nz + kk] * h[idx3(ii, jj + 1, kk, g.ny, g.nz)];
      if (kk > 0) v -= Tz[ii * (g.ny * (g.nz - 1)) + jj * (g.nz - 1) + (kk - 1)] * h[idx3(ii, jj, kk - 1, g.ny, g.nz)];
      if (kk < g.nz - 1) v -= Tz[ii * (g.ny * (g.nz - 1)) + jj * (g.nz - 1) + kk] * h[idx3(ii, jj, kk + 1, g.ny, g.nz)];
      out[c] = v + (alpha > 0 ? sum : 0);
    }
  };

  diagPC.set(diag);
  const x = pcg(applyOp, diagPC, b, { tol: 1e-11, maxIter: 20000 });

  const head = new Float64Array(x);
  if (!hasDirichlet) {
    const m = mean(head);
    for (let i = 0; i < n; i++) head[i] -= m;
  }

  // Расходы по граням (м3/сут) и скорости в ячейках
  const faceX = new Float64Array((g.nx - 1) * g.ny * g.nz);
  const faceY = new Float64Array(g.nx * (g.ny - 1) * g.nz);
  const faceZ = new Float64Array(g.nx * g.ny * (g.nz - 1));
  const darcyX = new Float64Array(n), darcyY = new Float64Array(n), darcyZ = new Float64Array(n);

  // Расходы считаются по гармоническому среднему проницаемости на грани —
  // той же величине, что используется при сборке матрицы.
  const kfx = facePerm(g, k, 0);
  const kfy = facePerm(g, k, 1);
  const kfz = facePerm(g, k, 2);
  for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx - 1; ii++) {
    const c = idx3(ii, jj, kk, g.ny, g.nz);
    const fi = ii * (g.ny * g.nz) + jj * g.nz + kk;
    const fl = -kfx[fi] * (head[idx3(ii + 1, jj, kk, g.ny, g.nz)] - head[c]) / g.dx;
    faceX[fi] = fl * g.dy * g.dz;
    darcyX[c] += 0.5 * fl;
    darcyX[idx3(ii + 1, jj, kk, g.ny, g.nz)] += 0.5 * fl;
  }
  for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny - 1; jj++) for (let ii = 0; ii < g.nx; ii++) {
    const c = idx3(ii, jj, kk, g.ny, g.nz);
    const fi = ii * ((g.ny - 1) * g.nz) + jj * g.nz + kk;
    const fl = -kfy[fi] * (head[idx3(ii, jj + 1, kk, g.ny, g.nz)] - head[c]) / g.dy;
    faceY[fi] = fl * g.dx * g.dz;
    darcyY[c] += 0.5 * fl;
    darcyY[idx3(ii, jj + 1, kk, g.ny, g.nz)] += 0.5 * fl;
  }
  for (let kk = 0; kk < g.nz - 1; kk++) for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx; ii++) {
    const c = idx3(ii, jj, kk, g.ny, g.nz);
    const fi = ii * (g.ny * (g.nz - 1)) + jj * (g.nz - 1) + kk;
    const fl = -kfz[fi] * kvRatio * (head[idx3(ii, jj, kk + 1, g.ny, g.nz)] - head[c]) / g.dz;
    faceZ[fi] = fl * g.dx * g.dy;
    darcyZ[c] += 0.5 * fl;
    darcyZ[idx3(ii, jj, kk + 1, g.ny, g.nz)] += 0.5 * fl;
  }

  const vx = new Float64Array(n), vy = new Float64Array(n), vz = new Float64Array(n);
  const speed = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const phi = porosity && porosity[i] > 1e-6 ? porosity[i] : 1e-6;
    vx[i] = darcyX[i] / phi; vy[i] = darcyY[i] / phi; vz[i] = darcyZ[i] / phi;
    speed[i] = Math.hypot(vx[i], vy[i], vz[i]);
  }

  // Баланс через границы
  let totalIn = 0, totalOut = 0, bIn = 0, bOut = 0;
  for (let i = 0; i < n; i++) {
    if (sources[i] > 0) totalIn += sources[i]; else totalOut -= sources[i];
  }
  for (const [face, cond] of Object.entries(bcFaces)) {
    if (cond !== 'dirichlet') continue;
    const [axis, side] = FACES[face];
    const hd = bcHeads[face] ?? 0;
    const half = axis === 0 ? 0.5 * g.dx : axis === 1 ? 0.5 * g.dy : 0.5 * g.dz;
    const area = axis === 0 ? g.dy * g.dz : axis === 1 ? g.dx * g.dz : g.dx * g.dy;
    for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx; ii++) {
      let onFace = false;
      if (axis === 0) onFace = (side === 0 ? ii === 0 : ii === g.nx - 1);
      else if (axis === 1) onFace = (side === 0 ? jj === 0 : jj === g.ny - 1);
      else onFace = (side === 0 ? kk === 0 : kk === g.nz - 1);
      if (!onFace) continue;
      const c = idx3(ii, jj, kk, g.ny, g.nz);
      const kEff = k[c] * (axis === 2 ? kvRatio : 1);
      const qIn = kEff * (hd - head[c]) / half * area;
      if (qIn > 0) bIn += qIn; else bOut += -qIn;
    }
  }

  return {
    head, vx, vy, vz, speed, darcyX, darcyY, darcyZ,
    faceX, faceY, faceZ,
    residual: totalIn - totalOut - (bIn - bOut),
    totalIn, totalOut, boundaryIn: bIn, boundaryOut: bOut,
    solver: 'pcg',
  };
}

/**
 * Сопряжённые градиенты с предобусловливанием Жаку.
 *
 * Для чистых условий Неймана матрица вырождена (константа — нулевое
 * пространство), поэтому применяется проекционная схема: итераты
 * удерживаются в подпространстве, ортогональном константе. Это позволяет
 * сойтись без потери точности, в отличие от простой регуляризации.
 */
function pcg(applyOp, diag, b, opts) {
  const n = b.length;
  const project = !!opts.project;
  const x = new Float64Array(n);
  const r = new Float64Array(b);
  const z = new Float64Array(n);
  const p = new Float64Array(n);
  const Ap = new Float64Array(n);

  const subtractMean = (v) => {
    if (!project) return;
    const m = mean(v);
    for (let i = 0; i < n; i++) v[i] -= m;
  };

  if (project) {
    subtractMean(r);
    // согласованность: среднее правой части должно быть нулём
    let m = mean(r);
    if (Math.abs(m) > 1e-12) for (let i = 0; i < n; i++) r[i] -= m;
  }

  let bnorm = 0;
  for (let i = 0; i < n; i++) bnorm += r[i] * r[i];
  bnorm = Math.sqrt(bnorm);
  if (bnorm < 1e-300) return x;

  for (let i = 0; i < n; i++) z[i] = r[i] / diag[i];
  subtractMean(z);
  for (let i = 0; i < n; i++) p[i] = z[i];
  let rz = 0;
  for (let i = 0; i < n; i++) rz += r[i] * z[i];
  const tol = opts.tol ?? 1e-11;
  let iter = 0;
  for (; iter < opts.maxIter; iter++) {
    applyOp(p, Ap);
    if (project) subtractMean(Ap);
    let pAp = 0;
    for (let i = 0; i < n; i++) pAp += p[i] * Ap[i];
    if (!(Math.abs(pAp) > 1e-300)) break;
    const alpha = rz / pAp;
    let rn = 0;
    for (let i = 0; i < n; i++) {
      x[i] += alpha * p[i];
      r[i] -= alpha * Ap[i];
    }
    if (project) subtractMean(r);
    for (let i = 0; i < n; i++) rn += r[i] * r[i];
    if (Math.sqrt(rn) / bnorm < tol) { iter++; break; }
    for (let i = 0; i < n; i++) z[i] = r[i] / diag[i];
    subtractMean(z);
    let rzNew = 0;
    for (let i = 0; i < n; i++) rzNew += r[i] * z[i];
    const beta = rzNew / rz;
    for (let i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
    if (project) subtractMean(p);
    rz = rzNew;
  }
  if (opts.iterations) opts.iterations(iter);
  return x;
}

// ------------------------------------------------------------------ линии тока

function trilinear(field, g, px, py, pz) {
  const fx = clamp((px - g.x0) / g.dx - 0.5, 0, g.nx - 1.000001);
  const fy = clamp((py - g.y0) / g.dy - 0.5, 0, g.ny - 1.000001);
  const fz = clamp((pz - g.z0) / g.dz - 0.5, 0, g.nz - 1.000001);
  const i0 = Math.floor(fx), j0 = Math.floor(fy), k0 = Math.floor(fz);
  const i1 = Math.min(i0 + 1, g.nx - 1), j1 = Math.min(j0 + 1, g.ny - 1), k1 = Math.min(k0 + 1, g.nz - 1);
  const tx = fx - i0, ty = fy - j0, tz = fz - k0;
  const at = (i, j, k) => field[idx3(i, j, k, g.ny, g.nz)];
  const c00 = at(i0, j0, k0) * (1 - tx) + at(i1, j0, k0) * tx;
  const c10 = at(i0, j1, k0) * (1 - tx) + at(i1, j1, k0) * tx;
  const c01 = at(i0, j0, k1) * (1 - tx) + at(i1, j0, k1) * tx;
  const c11 = at(i0, j1, k1) * (1 - tx) + at(i1, j1, k1) * tx;
  const c0 = c00 * (1 - ty) + c10 * ty;
  const c1 = c01 * (1 - ty) + c11 * ty;
  return c0 * (1 - tz) + c1 * tz;
}

export function traceStreamlines(g, vx, vy, vz, seeds, labels, stopCells, opts = {}) {
  const step = 0.35 * Math.min(g.dx, g.dy, g.dz);
  const stagnation = opts.stagnationThreshold ?? 0.05;
  const maxResidence = opts.maxResidence ?? 3650;
  const maxSteps = opts.maxSteps ?? 4000;
  const N = seeds.length;
  if (!N) return [];

  const P = seeds.map((s) => s.slice());
  const alive = new Uint8Array(N).fill(1);
  const outWell = new Array(N).fill('');
  const stagnant = new Uint8Array(N);
  const paths = Array.from({ length: N }, () => []);
  const speeds = Array.from({ length: N }, () => []);
  const times = Array.from({ length: N }, () => []);
  const curT = new Float64Array(N);

  for (let i = 0; i < N; i++) {
    paths[i].push(P[i].slice());
    speeds[i].push(Math.hypot(trilinear(vx, g, ...P[i]), trilinear(vy, g, ...P[i]), trilinear(vz, g, ...P[i])));
    times[i].push(0);
  }

  for (let step_i = 0; step_i < maxSteps; step_i++) {
    let anyAlive = false;
    for (let i = 0; i < N; i++) {
      if (!alive[i]) continue;
      anyAlive = true;
      const p = P[i];
      const vel = (q) => {
        const a = trilinear(vx, g, ...q), b = trilinear(vy, g, ...q), c = trilinear(vz, g, ...q);
        const m = Math.hypot(a, b, c);
        return m < 1e-14 ? [0, 0, 0] : [a / m, b / m, c / m];
      };
      const k1 = vel(p);
      const mid1 = [p[0] + 0.5 * step * k1[0], p[1] + 0.5 * step * k1[1], p[2] + 0.5 * step * k1[2]];
      const k2 = vel(mid1);
      const mid2 = [p[0] + 0.5 * step * k2[0], p[1] + 0.5 * step * k2[1], p[2] + 0.5 * step * k2[2]];
      const k3 = vel(mid2);
      const end = [p[0] + step * k3[0], p[1] + step * k3[1], p[2] + step * k3[2]];
      const k4 = vel(end);
      const np = [
        p[0] + step / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
        p[1] + step / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]),
        p[2] + step / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]),
      ];
      const sp = Math.hypot(trilinear(vx, g, ...np), trilinear(vy, g, ...np), trilinear(vz, g, ...np));
      const speedNow = Math.hypot(trilinear(vx, g, ...p), trilinear(vy, g, ...p), trilinear(vz, g, ...p));
      P[i] = np;
      curT[i] += step / Math.max(speedNow, 1e-9);
      paths[i].push(np);
      speeds[i].push(sp);
      times[i].push(curT[i]);

      const out = np[0] < g.x0 - 0.5 * g.dx || np[0] > g.x0 + g.nx * g.dx + 0.5 * g.dx
        || np[1] < g.y0 - 0.5 * g.dy || np[1] > g.y0 + g.ny * g.dy + 0.5 * g.dy
        || np[2] < g.z0 - 0.5 * g.dz || np[2] > g.z0 + g.nz * g.dz + 0.5 * g.dz;
      if (out) { alive[i] = 0; continue; }
      if (sp < stagnation) { stagnant[i] = 1; alive[i] = 0; continue; }
      if (curT[i] > maxResidence) { alive[i] = 0; continue; }
      const ci = clamp(Math.round((np[0] - g.x0) / g.dx - 0.5), 0, g.nx - 1);
      const cj = clamp(Math.round((np[1] - g.y0) / g.dy - 0.5), 0, g.ny - 1);
      const ck = clamp(Math.round((np[2] - g.z0) / g.dz - 0.5), 0, g.nz - 1);
      const key = idx3(ci, cj, ck, g.ny, g.nz);
      if (stopCells.has(key)) { outWell[i] = stopCells.get(key); alive[i] = 0; }
    }
    if (!anyAlive) break;
  }

  const result = [];
  for (let i = 0; i < N; i++) {
    if (paths[i].length < 2) continue;
    result.push({
      points: paths[i], speeds: speeds[i], times: times[i],
      wellIn: labels[i] ?? '', wellOut: outWell[i],
      stagnant: !!stagnant[i], residence: curT[i],
    });
  }
  return result;
}

export function seedsFromWells(wells, g) {
  const pts = [], labels = [];
  for (const w of wells) {
    if (w.well_type !== 'injection') continue;
    for (const c of wellScreenCells(w, g)) {
      // обратное преобразование плоского индекса (порядок C)
      const k = c % g.nz;
      const j = Math.floor(c / g.nz) % g.ny;
      const i = Math.floor(c / (g.ny * g.nz));
      pts.push([g.x0 + (i + 0.5) * g.dx, g.y0 + (j + 0.5) * g.dy, g.z0 + (k + 0.5) * g.dz]);
      labels.push(w.id);
    }
  }
  return { pts, labels };
}

export function stopCellsFromWells(wells, g) {
  const map = new Map();
  for (const w of wells) {
    if (w.well_type !== 'production') continue;
    for (const c of wellScreenCells(w, g)) map.set(c, w.id);
  }
  return map;
}

// ------------------------------------------------------------------ кинетика

export function temperatureFactor(params) {
  const t = (params.temperature ?? 20) + 273.15;
  const ea = (params.activation_energy ?? 20) * 1000;
  if (Math.abs(t - T0_REF) < 1e-6) return 1;
  return Math.exp(-ea / R_GAS * (1 / t - 1 / T0_REF));
}

function reactionRate(cAcid, cFe3, uRel, feoh3, params, reactions) {
  const tf = temperatureFactor(params);
  const colmat = 1 / (1 + (params.beta ?? 0) * Math.max(feoh3, 0));
  const solid = Math.cbrt(Math.max(uRel, 0) + 1e-12);
  let r = 0;
  if (reactions[0]) r += (params.k1 ?? 0) * Math.max(cAcid, 0) * solid;
  if (reactions[1]) r += (params.k2 ?? 0) * Math.max(cFe3, 0) * solid;
  return r * colmat * tf;
}

/**
 * Массоперенос в 3D: адвекция (upwind) + диффузия-дисперсия + реакции.
 * Полностью соответствует backend/app/chemistry.py.
 */
export function simulate3D(g, kArr, porosity, uranium, clay, faceX, faceY, faceZ,
  wells, params, outputTimes, reactions = [true, true, true, true], cfl = 0.6) {
  const n = g.n;
  const phi = new Float64Array(n);
  for (let i = 0; i < n; i++) phi[i] = porosity[i] > 1e-6 ? porosity[i] : 1e-6;
  const poreVol = new Float64Array(n);
  for (let i = 0; i < n; i++) poreVol[i] = phi[i] * g.cellVolume;

  const cAcid = new Float64Array(n).fill(params.acid_concentration ?? 20);
  const cU = new Float64Array(n);
  const cFe3 = new Float64Array(n).fill(params.fe3_conc ?? 0);
  const oreFrac = new Float64Array(n);
  for (let i = 0; i < n; i++) oreFrac[i] = Math.max((uranium[i] ?? 0) / 1000, 0);
  const uRel = new Float64Array(n).fill(1);
  const uRel0 = new Float64Array(n).fill(1);

  // Ячейки закачки и отбора с удельным расходом
  const inj = new Map(), prod = new Map();
  for (const w of wells) {
    const cells = wellScreenCells(w, g);
    if (!cells.length) continue;
    if (w.well_type === 'injection') {
      const per = w.rate / cells.length;
      for (const c of cells) inj.set(c, { acid: w.acid_in ?? 20, well: w.id, rate: per });
    } else if (w.well_type === 'production') {
      const per = Math.abs(w.rate) / cells.length;
      for (const c of cells) prod.set(c, { well: w.id, rate: per });
    }
  }
  for (const [c, info] of inj) cAcid[c] = info.acid;

  const totalTime = params.integration_days ?? 365;
  let outTimes = Array.from(new Set([...(outputTimes ?? []), totalTime])).sort((a, b) => a - b)
    .filter((t) => t <= totalTime + 1e-9);
  if (!outTimes.length) outTimes = [totalTime];
  if (outTimes[outTimes.length - 1] < totalTime - 1e-9) outTimes.push(totalTime);

  const res = {
    times: [], acid: [], uranium: [], residual: [], urRel: [],
    prUranium: [], vrAcid: [], cumulativeU: [], cumulativeAcid: [],
    prodVolume: [], injVolume: [],
  };

  const oreMask = new Uint8Array(n);
  for (let i = 0; i < n; i++) oreMask[i] = oreFrac[i] > 0 ? 1 : 0;

  let cumU = 0, cumAcid = 0, prodVol = 0, injVol = 0, prevU = 0, prevVol = 0;

  const record = (t) => {
    const dU = cumU - prevU, dV = prodVol - prevVol;
    const prU = dV > 1e-9 ? dU * 1000 / dV : 0;
    let vrA = 0, cnt = 0;
    for (let i = 0; i < n; i++) if (oreMask[i]) { vrA += cAcid[i]; cnt++; }
    vrA = cnt ? vrA / cnt : mean(cAcid);
    res.times.push(t);
    res.acid.push(cAcid.slice());
    res.uranium.push(cU.slice());
    const resArr = new Float64Array(n);
    for (let i = 0; i < n; i++) resArr[i] = uranium[i] * uRel[i];
    res.residual.push(resArr);
    res.urRel.push(uRel.slice());
    res.prUranium.push(prU);
    res.vrAcid.push(vrA);
    res.cumulativeU.push(cumU);
    res.cumulativeAcid.push(cumAcid);
    res.prodVolume.push(prodVol);
    res.injVolume.push(injVol);
    prevU = cumU; prevVol = prodVol;
  };
  record(0);

  // --- устойчивость ---
  const alphaL = 0.05, dMol = 0.01;
  const diffX = new Float64Array(faceX.length);
  const diffY = new Float64Array(faceY.length);
  const diffZ = new Float64Array(faceZ.length);
  const areaX = g.dy * g.dz, areaY = g.dx * g.dz, areaZ = g.dx * g.dy;
  for (let i = 0; i < faceX.length; i++) diffX[i] = (alphaL * Math.abs(faceX[i]) / areaX + dMol) * areaX / g.dx;
  for (let i = 0; i < faceY.length; i++) diffY[i] = (alphaL * Math.abs(faceY[i]) / areaY + dMol) * areaY / g.dy;
  for (let i = 0; i < faceZ.length; i++) diffZ[i] = (alphaL * Math.abs(faceZ[i]) / areaZ + dMol) * areaZ / g.dz;

  const stab = new Float64Array(n);
  const addStab = (cell, v) => { stab[cell] += v; };
  for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx; ii++) {
    const c = idx3(ii, jj, kk, g.ny, g.nz);
    if (ii > 0) { const f = faceX[(ii - 1) * (g.ny * g.nz) + jj * g.nz + kk]; addStab(c, Math.max(f, 0)); }
    if (ii < g.nx - 1) { const f = faceX[ii * (g.ny * g.nz) + jj * g.nz + kk]; addStab(c, Math.max(-f, 0)); }
    if (jj > 0) { const f = faceY[ii * ((g.ny - 1) * g.nz) + (jj - 1) * g.nz + kk]; addStab(c, Math.max(f, 0)); }
    if (jj < g.ny - 1) { const f = faceY[ii * ((g.ny - 1) * g.nz) + jj * g.nz + kk]; addStab(c, Math.max(-f, 0)); }
    if (kk > 0) { const f = faceZ[ii * (g.ny * (g.nz - 1)) + jj * (g.nz - 1) + (kk - 1)]; addStab(c, Math.max(f, 0)); }
    if (kk < g.nz - 1) { const f = faceZ[ii * (g.ny * (g.nz - 1)) + jj * (g.nz - 1) + kk]; addStab(c, Math.max(-f, 0)); }
  }
  for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx; ii++) {
    const c = idx3(ii, jj, kk, g.ny, g.nz);
    if (ii > 0) addStab(c, 2 * diffX[(ii - 1) + (g.nx - 1) * (jj + g.ny * kk)]);
    if (ii < g.nx - 1) addStab(c, 2 * diffX[ii + (g.nx - 1) * (jj + g.ny * kk)]);
    if (jj > 0) addStab(c, 2 * diffY[ii + g.nx * ((jj - 1) + (g.ny - 1) * kk)]);
    if (jj < g.ny - 1) addStab(c, 2 * diffY[ii * ((g.ny - 1) * g.nz) + jj * g.nz + kk]);
    if (kk > 0) addStab(c, 2 * diffZ[ii * (g.ny * (g.nz - 1)) + jj * (g.nz - 1) + (kk - 1)]);
    if (kk < g.nz - 1) addStab(c, 2 * diffZ[ii * (g.ny * (g.nz - 1)) + jj * (g.nz - 1) + kk]);
  }

  const dAcid = new Float64Array(n), dU = new Float64Array(n), dFe = new Float64Array(n);

  const flux = (F, cLo, cHi, dcoef) => {
    let f = 0;
    if (F > 0) f = F * cLo; else f = F * cHi;
    return f - dcoef * (cHi - cLo);
  };

  let t = 0, nextOut = 0, steps = 0;
  while (t < totalTime - 1e-9 && steps < 400000) {
    let dt = Infinity;
    for (let i = 0; i < n; i++) {
      if (stab[i] > 1e-12) dt = Math.min(dt, cfl * poreVol[i] / stab[i]);
    }
    if (!isFinite(dt)) dt = 1;
    dt = Math.min(dt, totalTime - t, 2);
    if (dt <= 0) break;

    dAcid.fill(0); dU.fill(0); dFe.fill(0);

    // Сбор потоков по граням
    for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx - 1; ii++) {
      const c0 = idx3(ii, jj, kk, g.ny, g.nz), c1 = idx3(ii + 1, jj, kk, g.ny, g.nz);
      const fi = ii * (g.ny * g.nz) + jj * g.nz + kk;
      const f = flux(faceX[fi], cAcid[c0], cAcid[c1], diffX[fi]);
      dAcid[c0] -= f * dt; dAcid[c1] += f * dt;
      const fu = flux(faceX[fi], cU[c0], cU[c1], diffX[fi]);
      dU[c0] -= fu * dt; dU[c1] += fu * dt;
      const ff = flux(faceX[fi], cFe3[c0], cFe3[c1], diffX[fi]);
      dFe[c0] -= ff * dt; dFe[c1] += ff * dt;
    }
    for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny - 1; jj++) for (let ii = 0; ii < g.nx; ii++) {
      const c0 = idx3(ii, jj, kk, g.ny, g.nz), c1 = idx3(ii, jj + 1, kk, g.ny, g.nz);
      const fi = ii * ((g.ny - 1) * g.nz) + jj * g.nz + kk;
      const f = flux(faceY[fi], cAcid[c0], cAcid[c1], diffY[fi]);
      dAcid[c0] -= f * dt; dAcid[c1] += f * dt;
      const fu = flux(faceY[fi], cU[c0], cU[c1], diffY[fi]);
      dU[c0] -= fu * dt; dU[c1] += fu * dt;
      const ff = flux(faceY[fi], cFe3[c0], cFe3[c1], diffY[fi]);
      dFe[c0] -= ff * dt; dFe[c1] += ff * dt;
    }
    for (let kk = 0; kk < g.nz - 1; kk++) for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx; ii++) {
      const c0 = idx3(ii, jj, kk, g.ny, g.nz), c1 = idx3(ii, jj, kk + 1, g.ny, g.nz);
      const fi = ii * (g.ny * (g.nz - 1)) + jj * (g.nz - 1) + kk;
      const f = flux(faceZ[fi], cAcid[c0], cAcid[c1], diffZ[fi]);
      dAcid[c0] -= f * dt; dAcid[c1] += f * dt;
      const fu = flux(faceZ[fi], cU[c0], cU[c1], diffZ[fi]);
      dU[c0] -= fu * dt; dU[c1] += fu * dt;
      const ff = flux(faceZ[fi], cFe3[c0], cFe3[c1], diffZ[fi]);
      dFe[c0] -= ff * dt; dFe[c1] += ff * dt;
    }

    // Реакции
    const feoh3 = params.feoh3_conc ?? 0;
    const uptakeArr = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let uptake = ORE_DENSITY * Math.max(1 - phi[i], 0) * oreFrac[i]
        * reactionRate(cAcid[i], cFe3[i], uRel[i], feoh3, params, reactions);
      if (uRel[i] <= 1e-8) uptake = 0;
      if (inj.has(i)) uptake = 0;
      const maxUp = uRel[i] * uRel0[i] * ORE_DENSITY * (1 - phi[i]) * oreFrac[i] / dt;
      uptake = Math.min(uptake, Math.max(maxUp, 0));
      uptakeArr[i] = uptake;
      dAcid[i] += (-ACID_PER_U * uptake) * dt * poreVol[i];
      dU[i] += (1000 * uptake) * dt * poreVol[i];
      dFe[i] += (-FE3_PER_U * uptake) * dt * poreVol[i];
    }

    // Потери урана (осаждение и сорбция)
    if (reactions[2] || reactions[3]) {
      for (let i = 0; i < n; i++) {
        let lossRate = 0;
        if (reactions[2]) lossRate += (params.k3 ?? 0) * cU[i];
        if (reactions[3]) lossRate += (params.k4 ?? 0) * cU[i] * clamp(clay[i], 0, 1);
        const lossMass = lossRate * dt * poreVol[i];
        dU[i] -= lossMass;
        const denom = Math.max(ORE_DENSITY * (1 - phi[i]) * oreFrac[i] * uRel0[i], 1e-12);
        const oreBack = (lossMass * 1e-3) / denom;
        if (uRel[i] > 0) uRel[i] = clamp(uRel[i] + oreBack, 0, 1);
      }
    }

    for (let i = 0; i < n; i++) {
      cAcid[i] = clamp(cAcid[i] + dAcid[i] / Math.max(poreVol[i], 1e-12), 0, 1e4);
      cU[i] = clamp(cU[i] + dU[i] / Math.max(poreVol[i], 1e-12), 0, 1e6);
      cFe3[i] = clamp(cFe3[i] + dFe[i] / Math.max(poreVol[i], 1e-12), 0, 1e4);
      // растворённый уран уменьшает остаточное содержание в руде
      const denom = Math.max(ORE_DENSITY * (1 - phi[i]) * oreFrac[i] * uRel0[i], 1e-12);
      uRel[i] = clamp(uRel[i] - uptakeArr[i] * dt / denom, 0, 1);
    }

    // расход кислоты на реакции, кг
    let acidUsedKg = 0;
    for (let i = 0; i < n; i++) acidUsedKg += ACID_PER_U * uptakeArr[i] * g.cellVolume * dt;

    // Отбор продуктивного раствора
    let stepVol = 0, stepMass = 0;
    for (const [c, info] of prod) {
      const q = info.rate;
      stepVol += q * dt;
      stepMass += cU[c] * 1e-3 * q * dt;
      const dMass = cU[c] * 1e-3 * q * dt;
      cU[c] = Math.max(0, cU[c] - dMass / Math.max(poreVol[c], 1e-12) * 1000);
      cAcid[c] = Math.max(0, cAcid[c] - (cAcid[c] * 1e-3 * q * dt) / Math.max(poreVol[c], 1e-12) * 1000);
      cFe3[c] = Math.max(0, cFe3[c] - (cFe3[c] * 1e-3 * q * dt) / Math.max(poreVol[c], 1e-12) * 1000);
    }
    prodVol += stepVol;
    cumU += stepMass;
    for (const [c, info] of inj) injVol += info.rate * dt;
    cumAcid += acidUsedKg;

    // Источники в закачивающих ячейках
    for (const [c, info] of inj) {
      cAcid[c] = info.acid;
      cU[c] = 0;
      cFe3[c] = params.fe3_conc ?? 0;
    }

    t += dt;
    steps++;
    if (nextOut < outTimes.length && t >= outTimes[nextOut] - 1e-9) {
      record(outTimes[nextOut]);
      nextOut++;
    }
  }
  if (res.times.length <= 1 || res.times[res.times.length - 1] < totalTime - 1e-6) record(totalTime);

  return res;
}

/** Кинетика вдоль линии тока (1D по времени пребывания). */
export function simulateFlowline(points, times, grade, porosity, params,
  reactions = [true, true, true, true]) {
  const n = points.length;
  if (n < 2) return { acid: [], uranium: [], rel: [], times: [] };
  let cAcid = params.flowline_acid ?? params.acid_concentration ?? 20;
  let cU = 0, cFe3 = params.fe3_conc ?? 0, uRel = 1;
  const acidP = [cAcid], uP = [cU], relP = [uRel], tP = [times[0]];
  const feoh3 = params.feoh3_conc ?? 0;

  for (let i = 0; i < n - 1; i++) {
    const dt = Math.max(0, times[i + 1] - times[i]);
    if (dt <= 0) { acidP.push(cAcid); uP.push(cU); relP.push(uRel); tP.push(times[i + 1]); continue; }
    const steps = Math.max(1, Math.min(200, Math.max(2, dt * 4)));
    const sdt = dt / steps;
    const oreFrac = Math.max((grade[i] ?? 0) / 1000, 0);
    const phi = porosity[i] ?? 0.15;
    for (let s = 0; s < steps; s++) {
      let uptake = ORE_DENSITY * Math.max(1 - phi, 0) * oreFrac
        * reactionRate(cAcid, cFe3, uRel, feoh3, params, reactions);
      const maxUp = uRel * ORE_DENSITY * (1 - phi) * oreFrac / sdt;
      uptake = Math.min(uptake, Math.max(maxUp, 0));
      cAcid = Math.max(0, cAcid + (-ACID_PER_U * uptake) * sdt);
      cU = Math.max(0, cU + (1000 * uptake) * sdt);
      cFe3 = Math.max(0, cFe3 + (-FE3_PER_U * uptake) * sdt);
      const denom = Math.max(ORE_DENSITY * (1 - phi) * oreFrac, 1e-12);
      uRel = Math.max(0, uRel - uptake * sdt / denom);
    }
    acidP.push(cAcid); uP.push(cU); relP.push(uRel); tP.push(times[i + 1]);
  }
  return { acid: acidP, uranium: uP, rel: relP, times: tP };
}

// ------------------------------------------------------------------ запасы

export function calculateReserves(g, grade, acid, opts = {}) {
  const threshold = opts.acidThreshold ?? 5;
  const contour = opts.contourMask ?? null;
  const density = opts.density ?? 2200;
  const n = g.n ?? g.nx * g.ny * g.nz;
  const cellOreT = density * g.cellVolume / 1000;
  const res = {
    reserveKg: 0, reserveT: 0, acidifiedVolume: 0, acidifiedCells: 0,
    areaM2: 0, areaHa: 0, avgGrade: 0, avgGradePct: 0, avgThickness: 0,
    weightedMetrePercent: 0, minGrade: 0, maxGrade: 0, totalOreT: 0,
    contourAreaM2: 0, contourVolume: 0, byLayer: [],
  };

  let inContour = 0;
  for (let i = 0; i < n; i++) if (!contour || contour[i]) inContour++;
  res.totalOreT = inContour * cellOreT;

  let wSum = 0, gSum = 0, minG = Infinity, maxG = -Infinity;
  for (let i = 0; i < n; i++) {
    if (contour && !contour[i]) continue;
    if (acid[i] < threshold) continue;
    res.acidifiedCells++;
    const gi = grade[i];
    wSum += cellOreT;
    gSum += gi * cellOreT;
    if (gi < minG) minG = gi;
    if (gi > maxG) maxG = gi;
  }
  res.acidifiedVolume = res.acidifiedCells * g.cellVolume;
  if (!res.acidifiedCells) { res.minGrade = 0; res.maxGrade = 0; return res; }
  res.reserveKg = gSum;
  res.reserveT = res.reserveKg / 1000;
  res.avgGrade = gSum / wSum;
  res.avgGradePct = res.avgGrade / 10;
  res.minGrade = minG;
  res.maxGrade = maxG;

  // Мощность и площадь в плане
  const activeZ = new Float64Array(g.nx * g.ny);
  const oreThk = new Float64Array(g.nx * g.ny);
  let cols = 0;
  for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx; ii++) {
    const i = idx3(ii, jj, kk, g.ny, g.nz);
    if (contour && !contour[i]) continue;
    if (acid[i] < threshold) continue;
    const c = ii * g.ny + jj;
    activeZ[c] += g.dz;
    oreThk[c] += grade[i];
  }
  for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx; ii++) {
    const c = ii * g.ny + jj;
    if (activeZ[c] > 0) { cols++; oreThk[c] /= (activeZ[c] / g.dz); }
  }
  res.avgThickness = activeZ.length ? activeZ.reduce((a, b) => a + b, 0) / activeZ.length : 0;
  let num = 0, den = 0;
  for (let c = 0; c < activeZ.length; c++) {
    if (activeZ[c] > 0) { num += oreThk[c] * activeZ[c]; den += activeZ[c]; }
  }
  res.weightedMetrePercent = den > 0 ? num / den : 0;
  res.areaM2 = cols * g.dx * g.dy;
  res.areaHa = res.areaM2 / 10000;

  for (let kk = 0; kk < g.nz; kk++) {
    let cnt = 0, gradeSum = 0, acidSum = 0;
    for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx; ii++) {
      const i = idx3(ii, jj, kk, g.ny, g.nz);
      if (contour && !contour[i]) continue;
      if (acid[i] < threshold) continue;
      cnt++; gradeSum += grade[i]; acidSum += acid[i];
    }
    if (!cnt) continue;
    res.byLayer.push({
      layer: kk + 1, zTop: g.z0 + kk * g.dz, zBottom: g.z0 + (kk + 1) * g.dz,
      cells: cnt, volume: cnt * g.cellVolume,
      reserveKg: (gradeSum / cnt) * cnt * cellOreT,
      avgGrade: gradeSum / cnt, avgAcid: acidSum / cnt,
    });
  }
  return res;
}

export function contourFromPolygon(g, polygon) {
  const mask = new Uint8Array(g.n);
  if (!polygon || polygon.length < 3) return mask;
  const c = gridCenters(g);
  for (let kk = 0; kk < g.nz; kk++) for (let jj = 0; jj < g.ny; jj++) for (let ii = 0; ii < g.nx; ii++) {
    if (pointInPolygon(c.x[ii], c.y[jj], polygon)) mask[idx3(ii, jj, kk, g.ny, g.nz)] = 1;
  }
  return mask;
}

function pointInPolygon(x, y, poly) {
  let inside = false;
  for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
    const xi = poly[i][0], yi = poly[i][1], xj = poly[j][0], yj = poly[j][1];
    if ((yi > y) !== (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi + 1e-30) + xi) inside = !inside;
  }
  return inside;
}

// ------------------------------------------------------------------ импорт

const FIELD_ALIASES = {
  name: ['скважина', 'наименование', 'название', 'well', 'wellname', 'name', '№ скважины', 'номер скважины', 'скв', 'скв.', 'id', 'borehole'],
  x: ['x', 'координата x', 'координатах', 'координата_х', 'восток', 'восточный', 'coord_x', 'coordinate_x', 'easting'],
  y: ['y', 'координата y', 'координатау', 'координата_у', 'север', 'северный', 'coord_y', 'coordinate_y', 'northing'],
  z: ['z', 'глубина', 'отметка', 'абс отметка', 'абсолютная отметка', 'отметка фильтра', 'elev', 'elevation', 'coord_z', 'глубина фильтра'],
  k: ['кф', 'коэффициент фильтрации', 'коэффициент_фильтрации', 'проницаемость', 'k', 'permeability', 'kh', 'кп'],
  porosity: ['пористость', 'porosity', 'poro', 'п'],
  uranium: ['уран', 'содержание урана', 'cу', 'cu', 'u', 'uranium', 'grade', 'запас', 'содержание', 'уран, кг/т', 'u_pct'],
  clay: ['глина', 'доля глины', 'глинистость', 'clay', 'глин', 'песчаник'],
  lithology: ['литология', 'lithology', 'lit', 'тип породы', 'порода', 'lith'],
  rate: ['дебит', 'расход', 'q', 'rate', 'производительность', 'приток', 'деб', 'q_нагн'],
  well_type: ['тип', 'назначение', 'type', 'well_type', 'категория', 'скважина назначения'],
  acid_in: ['кислотность', 'кислота', 'h2so4', 'acid', 'вр', 'концентрация кислоты'],
  screen_top: ['верхний интервал', 'топ', 'top', 'начало интервала', 'screen_top', 'интервал перфорации верх', 'верх фильтра', 'верхний фильтр', 'верх отверстия', 'начало фильтра', 'верх перфорации', 'отметка верха фильтра'],
  screen_bottom: ['нижний интервал', 'низ', 'bottom', 'конец интервала', 'screen_bottom', 'интервал перфорации низ', 'низ фильтра', 'нижний фильтр', 'низ отверстия', 'конец фильтра', 'низ перфорации', 'отметка низа фильтра'],
};
const NUMERIC_FIELDS = new Set(['x', 'y', 'z', 'k', 'porosity', 'uranium', 'clay', 'rate', 'acid_in', 'screen_top', 'screen_bottom']);

const norm = (s) => String(s).trim().toLowerCase().replace(/[^a-zа-я0-9]+/g, '_').replace(/^_|_$/g, '');

export function autoMapColumns(columns) {
  const normCols = new Map(columns.map((c) => [norm(c), c]));
  const mapping = {};
  const used = new Set();
  for (const [field, aliases] of Object.entries(FIELD_ALIASES)) {
    let found = null;
    for (const alias of aliases) {
      const a = norm(alias);
      if (normCols.has(a) && !used.has(normCols.get(a))) { found = normCols.get(a); break; }
    }
    if (!found) {
      for (const [nc, orig] of normCols) {
        if (used.has(orig)) continue;
        for (const alias of aliases) {
          const a = norm(alias);
          if (a && (nc.includes(a) || a.includes(nc)) && a.length >= 3) { found = orig; break; }
        }
        if (found) break;
      }
    }
    mapping[field] = found ?? null;
    if (found) used.add(found);
  }
  return mapping;
}

function toFloat(v) {
  if (v == null || v === '') return NaN;
  if (typeof v === 'number') return v;
  const s = String(v).trim().replace(',', '.');
  const m = s.match(/[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?/);
  return m ? parseFloat(m[0]) : NaN;
}

function normType(v) {
  const s = String(v ?? '').trim().toLowerCase();
  if (['закачка', 'закачивающая', 'нагн', 'нагнетательная', 'injection', 'зак'].includes(s)) return 'injection';
  if (['откачка', 'откачивающая', 'добыча', 'отбор', 'production', 'отк'].includes(s)) return 'production';
  if (['наблюдательная', 'контрольная', 'observation', 'набл'].includes(s)) return 'observation';
  return 'injection';
}

/** Разбор CSV-текста с автоопределением разделителя. */
export function parseDelimited(text) {
  const clean = text.replace(/^\uFEFF/, '');
  const firstLine = clean.split(/\r?\n/)[0] ?? '';
  const counts = {
    ';': (firstLine.match(/;/g) || []).length,
    ',': (firstLine.match(/,/g) || []).length,
    '\t': (firstLine.match(/\t/g) || []).length,
  };
  let sep = ';';
  let best = counts[';'];
  for (const c of [',', '\t']) if (counts[c] > best) { best = counts[c]; sep = c; }
  const rows = [];
  let row = [], field = '', inQuotes = false;
  for (let i = 0; i < clean.length; i++) {
    const ch = clean[i];
    if (inQuotes) {
      if (ch === '"') {
        if (clean[i + 1] === '"') { field += '"'; i++; } else inQuotes = false;
      } else field += ch;
    } else if (ch === '"') inQuotes = true;
    else if (ch === sep) { row.push(field); field = ''; }
    else if (ch === '\n') { row.push(field); rows.push(row); row = []; field = ''; }
    else if (ch !== '\r') field += ch;
  }
  if (field.length || row.length) { row.push(field); rows.push(row); }
  if (!rows.length) return { columns: [], data: [] };
  const columns = rows[0].map((c) => c.trim());
  const data = rows.slice(1).filter((r) => r.some((v) => String(v).trim() !== '')).map((r) => {
    const o = {};
    columns.forEach((c, i) => { o[c] = r[i] ?? ''; });
    return o;
  });
  return { columns, data };
}

export function rowsToWells(data, mapping) {
  const out = [];
  for (const row of data) {
    const rec = {};
    for (const [field, col] of Object.entries(mapping)) {
      if (!col) { rec[field] = null; continue; }
      const raw = row[col];
      if (NUMERIC_FIELDS.has(field)) rec[field] = toFloat(raw);
      else if (field === 'lithology') rec[field] = raw == null ? null : String(raw).trim();
      else if (field === 'well_type') rec[field] = normType(raw);
      else rec[field] = raw == null ? null : String(raw).trim();
    }
    if (!Number.isFinite(rec.x) || !Number.isFinite(rec.y)) continue;
    out.push(rec);
  }
  return out;
}
