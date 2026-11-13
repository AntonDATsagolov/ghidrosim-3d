/**
 * Сверка JS-движка с Python-эталоном.
 * Запуск:  node scripts/validate_engine.mjs <input.json> <output.json>
 *
 * input.json — проект (grid, wells, props, chem, bc) и опции расчёта;
 * output.json — результаты для сравнения с эталоном.
 */
import { readFileSync, writeFileSync } from 'node:fs';
import {
  makeGrid, buildPropertyModels, buildSources, solveDarcy, traceStreamlines,
  seedsFromWells, stopCellsFromWells, simulate3D, calculateReserves,
} from '../frontend/js/engine.js';

const inputPath = process.argv[2];
const outputPath = process.argv[3];
if (!inputPath || !outputPath) {
  console.error('Использование: node validate_engine.mjs <input.json> <output.json>');
  process.exit(2);
}

const input = JSON.parse(readFileSync(inputPath, 'utf8'));
const project = input.project;
const opts = input.options ?? {};

const t0 = Date.now();
const grid = makeGrid(project.grid);
const props = buildPropertyModels(project, opts);

const { q: sources } = buildSources(grid, project.wells);
const hydro = solveDarcy(
  grid, props.k, sources, props.porosity,
  project.bcFaces ?? {}, project.bcHeads ?? {}, opts.kv_ratio ?? 1.0);

const { pts, labels } = seedsFromWells(project.wells, grid);
const stops = stopCellsFromWells(project.wells, grid);
const lines = traceStreamlines(grid, hydro.vx, hydro.vy, hydro.vz, pts, labels, stops, {
  stagnationThreshold: opts.stagnation_threshold ?? 0.05,
  maxResidence: Math.max(365, (project.chem.integration_days ?? 365) * 3),
});

const times = [];
for (let t = opts.output_every ?? 30; t <= (project.chem.integration_days ?? 365); t += opts.output_every ?? 30) times.push(t);
const mt = simulate3D(grid, props.k, props.porosity, props.uranium, props.clay,
  hydro.faceX, hydro.faceY, hydro.faceZ, project.wells, project.chem, times,
  opts.reactions ?? [true, true, true, true]);

const reserves = calculateReserves(grid, mt.residual[mt.residual.length - 1],
  mt.acid[mt.acid.length - 1], { acidThreshold: opts.acid_threshold ?? 5 });

// Сводка для сравнения
const out = {
  elapsed_ms: Date.now() - t0,
  hydro: {
    residual: hydro.residual,
    total_in: hydro.totalIn,
    total_out: hydro.totalOut,
    head_min: Math.min(...hydro.head),
    head_max: Math.max(...hydro.head),
    speed_max: Math.max(...hydro.speed),
    speed_mean: hydro.speed.reduce((a, b) => a + b, 0) / hydro.speed.length,
  },
  property: {
    k_min: Math.min(...props.k), k_max: Math.max(...props.k),
    k_mean: props.k.reduce((a, b) => a + b, 0) / props.k.length,
    u_mean: props.uranium.reduce((a, b) => a + b, 0) / props.uranium.length,
  },
  streamlines: {
    count: lines.length,
    reached: lines.filter((l) => l.wellOut).length,
    residence_median: median(lines.map((l) => l.residence)),
  },
  chemistry: {
    times: mt.times,
    pr_uranium: mt.prUranium,
    vr_acid: mt.vrAcid,
    cumulative_u: mt.cumulativeU,
    cumulative_acid: mt.cumulativeAcid,
    prod_volume: mt.prodVolume,
    acid_final_mean: mean(mt.acid[mt.acid.length - 1]),
    uranium_final_max: Math.max(...mt.uranium[mt.uranium.length - 1]),
    residual_mean: mean(mt.residual[mt.residual.length - 1]),
  },
  reserves: {
    reserve_kg: reserves.reserveKg,
    reserve_t: reserves.reserveT,
    acidified_cells: reserves.acidifiedCells,
    acidified_volume: reserves.acidifiedVolume,
    area_ha: reserves.areaHa,
    avg_grade: reserves.avgGrade,
    avg_thickness: reserves.avgThickness,
    weighted_metre_percent: reserves.weightedMetrePercent,
  },
};

writeFileSync(outputPath, JSON.stringify(out, null, 1));
console.log(`JS-расчёт завершён за ${out.elapsed_ms} мс; невязка ${out.hydro.residual.toExponential(2)}`);

function mean(a) { return a.reduce((x, y) => x + y, 0) / a.length; }
function median(a) {
  if (!a.length) return 0;
  const s = [...a].sort((x, y) => x - y);
  const m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}
