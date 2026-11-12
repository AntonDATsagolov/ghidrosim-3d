/** Графики на canvas без внешних зависимостей. */

const AXIS = '#3a4654';
const GRID = '#232c36';
const TEXT = '#a9b6c4';

function setup(canvas) {
  const dpr = Math.min(devicePixelRatio || 1, 2);
  const w = canvas.clientWidth || canvas.width;
  const h = canvas.clientHeight || canvas.height;
  canvas.width = w * dpr;
  canvas.height = h * dpr;
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  return { ctx, w, h };
}

function niceTicks(min, max, count = 5) {
  if (!isFinite(min) || !isFinite(max)) return [0, 1];
  if (min === max) return [min];
  const span = max - min;
  const step0 = span / count;
  const mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const norm = step0 / mag;
  let step;
  if (norm <= 1) step = 1; else if (norm <= 2) step = 2;
  else if (norm <= 5) step = 5; else step = 10;
  step *= mag;
  const start = Math.ceil(min / step) * step;
  const out = [];
  for (let v = start; v <= max + step * 1e-6; v += step) out.push(v);
  return out;
}

/**
 * Линейный график с одной или двумя осями Y.
 * series: [{name, color, values, axis: 'left'|'right', dashed}]
 */
export function lineChart(canvas, series, opts = {}) {
  const { ctx, w, h } = setup(canvas);
  const padL = 58, padR = opts.dualAxis ? 52 : 16, padT = 14, padB = 34;
  const plotW = w - padL - padR;
  const plotH = h - padT - padB;
  if (plotW <= 10 || plotH <= 10) return;

  const xs = opts.x || (series[0] ? series[0].x || series[0].values.map((_, i) => i) : []);
  const xMin = opts.xMin ?? Math.min(...xs);
  const xMax = opts.xMax ?? Math.max(...xs);
  const left = series.filter((s) => s.axis !== 'right');
  const right = series.filter((s) => s.axis === 'right');
  const extent = (arr) => {
    let mn = Infinity, mx = -Infinity;
    arr.forEach((s) => s.values.forEach((v) => {
      if (v == null || !isFinite(v)) return;
      mn = Math.min(mn, v); mx = Math.max(mx, v);
    }));
    if (!isFinite(mn)) { mn = 0; mx = 1; }
    if (mn === mx) { mx = mn + 1; }
    return [mn, mx];
  };
  const [lMin, lMax] = left.length ? extent(left) : [0, 1];
  const [rMin, rMax] = right.length ? extent(right) : [0, 1];
  const lPad = (lMax - lMin) * 0.08;
  const rPad = (rMax - rMin) * 0.08;
  const L0 = lMin - lPad, L1 = lMax + lPad;
  const R0 = rMin - rPad, R1 = rMax + rPad;

  const px = (x) => padL + ((x - xMin) / (xMax - xMin || 1)) * plotW;
  const pyL = (v) => padT + plotH - ((v - L0) / (L1 - L0 || 1)) * plotH;
  const pyR = (v) => padT + plotH - ((v - R0) / (R1 - R0 || 1)) * plotH;

  // Сетка и подписи оси Y
  ctx.font = '11px ui-monospace, Menlo, monospace';
  ctx.strokeStyle = GRID;
  ctx.fillStyle = TEXT;
  ctx.lineWidth = 1;
  const ticks = niceTicks(L0, L1, 5);
  ticks.forEach((t) => {
    const y = pyL(t);
    ctx.beginPath();
    ctx.moveTo(padL, y);
    ctx.lineTo(padL + plotW, y);
    ctx.stroke();
    ctx.textAlign = 'right';
    ctx.fillText(fmt(t), padL - 6, y + 4);
  });
  if (right.length) {
    niceTicks(R0, R1, 5).forEach((t) => {
      const y = pyR(t);
      ctx.textAlign = 'left';
      ctx.fillStyle = '#7f8c9b';
      ctx.fillText(fmt(t), padL + plotW + 6, y + 4);
    });
  }

  // Ось X
  ctx.strokeStyle = AXIS;
  ctx.beginPath();
  ctx.moveTo(padL, padT);
  ctx.lineTo(padL, padT + plotH);
  ctx.lineTo(padL + plotW, padT + plotH);
  ctx.stroke();
  ctx.fillStyle = TEXT;
  ctx.textAlign = 'center';
  niceTicks(xMin, xMax, 6).forEach((t) => {
    const x = px(t);
    ctx.fillText(fmt(t), x, padT + plotH + 15);
  });
  if (opts.xLabel) {
    ctx.fillText(opts.xLabel, padL + plotW / 2, h - 6);
  }

  series.forEach((s) => {
    const py = s.axis === 'right' ? pyR : pyL;
    ctx.strokeStyle = s.color;
    ctx.lineWidth = s.width || 1.6;
    ctx.setLineDash(s.dashed ? [5, 4] : []);
    ctx.beginPath();
    let started = false;
    s.values.forEach((v, i) => {
      if (v == null || !isFinite(v)) { started = false; return; }
      const x = px(xs[i] ?? i);
      const y = py(v);
      if (!started) { ctx.moveTo(x, y); started = true; } else ctx.lineTo(x, y);
    });
    ctx.stroke();
    ctx.setLineDash([]);
  });
}

function fmt(v) {
  const a = Math.abs(v);
  if (a >= 1e6) return (v / 1e6).toFixed(1) + 'M';
  if (a >= 1e4) return (v / 1e3).toFixed(1) + 'k';
  if (a >= 100) return v.toFixed(0);
  if (a >= 1) return v.toFixed(1);
  if (a === 0) return '0';
  return v.toFixed(3);
}

/** Тепловая карта (план) с цветовой шкалой. */
export function heatmap(canvas, data, nx, ny, opts = {}) {
  const { ctx, w, h } = setup(canvas);
  const padL = 46, padB = 30, padT = 10, padR = 12;
  const plotW = w - padL - padR;
  const plotH = h - padT - padB;
  if (plotW <= 10 || plotH <= 10) return;
  const cw = plotW / nx;
  const ch = plotH / ny;
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const v = data[i * ny + j];
      const t = opts.norm ? opts.norm(v) : 0.5;
      ctx.fillStyle = opts.color ? opts.color(t) : defaultColor(t);
      ctx.fillRect(padL + i * cw, padT + (ny - 1 - j) * ch, cw + 0.6, ch + 0.6);
    }
  }
  ctx.strokeStyle = AXIS;
  ctx.strokeRect(padL, padT, plotW, plotH);
  ctx.fillStyle = TEXT;
  ctx.font = '11px ui-monospace, Menlo, monospace';
  ctx.textAlign = 'right';
  ctx.fillText(opts.yMax ?? '', padL - 6, padT + 10);
  ctx.fillText(opts.yMin ?? '', padL - 6, padT + plotH);
  ctx.textAlign = 'center';
  ctx.fillText(opts.xMin ?? '', padL, padT + plotH + 14);
  ctx.fillText(opts.xMax ?? '', padL + plotW, padT + plotH + 14);
  ctx.textAlign = 'left';
  ctx.fillText(opts.title || '', padL, padT - 1);
}

function defaultColor(t) {
  const stops = [
    [0.00, 27, 20, 100], [0.18, 31, 111, 178], [0.40, 35, 160, 138],
    [0.58, 90, 199, 106], [0.76, 217, 198, 46], [0.90, 232, 118, 44],
    [1.00, 200, 30, 30],
  ];
  t = Math.max(0, Math.min(1, t));
  const n = stops.length - 1;
  const i = Math.min(Math.floor(t * n), n - 1);
  const f = t * n - i;
  const a = stops[i], b = stops[i + 1];
  return `rgb(${Math.round(a[1] + (b[1] - a[1]) * f)},${Math.round(a[2] + (b[2] - a[2]) * f)},${Math.round(a[3] + (b[3] - a[3]) * f)})`;
}

/** Отрисовка легенды цветовой шкалы. */
export function drawLegend(canvas, min, max, colorFn) {
  const { ctx, w, h } = setup(canvas);
  for (let x = 0; x < w; x++) {
    const t = x / (w - 1);
    ctx.fillStyle = colorFn(t);
    ctx.fillRect(x, 0, 1, h);
  }
}
