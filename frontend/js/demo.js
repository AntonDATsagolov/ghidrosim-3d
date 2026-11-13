/**
 * Генератор демонстрационного проекта для автономной HTML-версии.
 * Данные синтетические: чередующиеся ряды закачных и откачных скважин,
 * две отработки, неоднородная проницаемость, рудные «линзы».
 */

/** Детерминированный ГПСЧ (mulberry32). */
function rng(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export function buildDemoProject() {
  const rand = rng(20240115);
  const grid = {
    x0: 0, y0: 0, z0: 0, nx: 32, ny: 24, nz: 6,
    dx: 25, dy: 25, dz: 10,
  };
  const wells = [];
  const nRows = 6, perRow = 7;
  const xStart = 100, xStep = 100, yStart = 100, yStep = 90;

  const zoneK = (x, y) => {
    let base = 1.6;
    if (x > 480) base *= 0.45;       // слабопроницаемый блок
    if (y > 250 && y < 330) base *= 1.8; // русло
    return base;
  };

  let idx = 0;
  for (let row = 0; row < nRows; row++) {
    const isInjection = row % 2 === 0;
    for (let c = 0; c < perRow; c++) {
      const x = xStart + c * xStep;
      const y = yStart + row * yStep;
      const zHead = 2.0 + (rand() - 0.5) * 3;
      idx += 1;
      const wid = 'W' + String(idx).padStart(2, '0');
      const screens = [{ top: 30, bottom: 45 }];
      if (row === 0 || row === nRows - 1) {
        screens.unshift({ top: 12, bottom: 24 });
      }
      const rate = (isInjection ? 620 : -600) * (0.9 + rand() * 0.2);
      wells.push({
        id: wid,
        name: (isInjection ? 'З' : 'О') + '-' + String(idx).padStart(2, '0'),
        well_type: isInjection ? 'injection' : 'production',
        trajectory: [{ x, y, z: zHead }],
        screens,
        screen: screens[0],
        rate,
        acid_in: isInjection ? 20 : 0,
        notes: 'Учебные данные',
      });
    }
  }

  // Баланс расходов: отбор равен закачке
  let injTotal = 0;
  for (const w of wells) if (w.well_type === 'injection') injTotal += w.rate;
  const prod = wells.filter((w) => w.well_type === 'production');
  if (prod.length) for (const w of prod) w.rate = -injTotal / prod.length;

  const permeability = [], uranium = [], porosity = [], clay = [], lithology = [];
  const push = (list, x, y, z, v) => list.push([x, y, z, v]);

  for (const w of wells) {
    const [x, y] = [w.trajectory[0].x, w.trajectory[0].y];
    for (const s of w.screens) {
      const zMid = 0.5 * (s.top + s.bottom);
      push(permeability, x, y, zMid, zoneK(x, y) * (0.75 + rand() * 0.55));
      let g = 0.55 * (0.4 + rand() * 1.2);
      if (x > 300 && x < 600 && y > 150 && y < 450) g *= 1.7;
      push(uranium, x, y, zMid, g);
      push(porosity, x, y, zMid, 0.13 + rand() * 0.13);
      push(clay, x, y, zMid, x > 480 ? 0.18 : 0.07);
      push(lithology, x, y, zMid, x > 480 ? 'глинистый алевролит' : 'песчаник');
    }
  }
  for (let i = 0; i < 40; i++) {
    const x = 30 + rand() * 740;
    const y = 30 + rand() * 540;
    const z = rand() < 0.5 ? 18 : 37;
    push(permeability, x, y, z, zoneK(x, y) * (0.8 + rand() * 0.45));
    push(uranium, x, y, z, 0.5 * (0.3 + rand() * 1.2));
    push(porosity, x, y, z, 0.13 + rand() * 0.12);
    push(clay, x, y, z, x > 480 ? 0.15 : 0.07);
    push(lithology, x, y, z, x > 480 ? 'глинистый алевролит' : 'песчаник');
  }

  return {
    name: 'Демо: блок 5-2-1 (учебные данные)',
    grid,
    wells,
    props: { permeability, porosity, uranium, clay, lithology },
    chem: {
      acid_concentration: 20, use_acid_data: true,
      fe3_conc: 2, feoh3_conc: 0.35, caco3_conc: 0,
      integration_days: 365,
      k1: 0.0003, k2: 0.00005, k3: 0.00002, k4: 0.00005,
      beta: 0.6, temperature: 22, activation_energy: 20,
      flowline_k1: 0.0003, flowline_k2: 0.00005, flowline_k3: 0.00002,
      flowline_k4: 0.00005, flowline_beta: 0.6, flowline_acid: 20,
    },
    bcFaces: {
      x_min: 'neumann', x_max: 'neumann', y_min: 'neumann',
      y_max: 'neumann', z_min: 'neumann', z_max: 'neumann',
    },
    bcHeads: {},
  };
}
