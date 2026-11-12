"""Модуль химической кинетики и массопереноса для ПСВ урана.

Модель процессов (упрощённая инженерная, калибруется по фактическим
кривым ПР/ВР):

    Основное реакция — сернокислотное выщелачивание уранинита:
        UO2 + H2SO4 -> UO2SO4 + H2O
    Реакция с Fe(III) — восстановление железа ураном (каталитический цикл):
        2Fe3+ + UO2 -> 2Fe2+ + U(VI)

Скорость растворения урана (смешанный контроль: реагент — 1-я порядок,
твердая фаза — 1/3):
        R = (k1 * C_a + k2 * C_Fe3) * U_rel^(1/3) / (1 + beta * C_Fe(OH)3)

Побочные потери урана:
    k3 — осаждение урана из раствора;
    k4 — сорбция урана на глинистой компоненте породы.

Температурная поправка — уравнение Аррениуса.

Стехиометрия (на 1 кг растворённого урана):
    a_acid  — 2.0 кг H2SO4 / кг U
    f_fe3   — 0.47 кг Fe(III) / кг U
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .models import ChemicalParams, GridSpec, WellType

R_GAS = 8.314  # Дж/(моль·К)
ORE_DENSITY = 2200.0  # кг/м3 — плотность руды с учётом скелета
ACID_PER_U = 2.0  # кг H2SO4 на кг U
FE3_PER_U = 0.47  # кг Fe(III) на кг U
T0_REF = 293.15  # К — температура, при которой заданы коэффициенты


@dataclass
class Kinetics:
    """Ядро кинетических уравнений (используется и в 3D, и в расчёте по линиям тока)."""

    params: ChemicalParams

    def temperature_factor(self) -> float:
        t = self.params.temperature + 273.15
        ea = self.params.activation_energy * 1000.0  # Дж/моль
        if abs(t - T0_REF) < 1e-6:
            return 1.0
        return float(np.exp(-ea / R_GAS * (1.0 / t - 1.0 / T0_REF)))

    def rate(self, c_acid: np.ndarray, c_fe3: np.ndarray, u_rel: np.ndarray,
             feoh3_pct: np.ndarray, reactions: tuple[bool, bool, bool, bool] = (True, True, True, True)) -> np.ndarray:
        """Суммарная скорость растворения урана, 1/сут."""
        p = self.params
        tf = self.temperature_factor()
        u_rel = np.clip(u_rel, 0.0, 1.0)
        colmat = 1.0 / (1.0 + p.beta * np.maximum(feoh3_pct, 0.0))
        solid_term = np.cbrt(np.maximum(u_rel, 0.0) + 1e-12)
        r = np.zeros_like(np.asarray(u_rel, dtype=float))
        if reactions[0]:
            r += p.k1 * np.maximum(c_acid, 0.0) * solid_term
        if reactions[1]:
            r += p.k2 * np.maximum(c_fe3, 0.0) * solid_term
        return r * colmat * tf

    def uranium_uptake_rate(self, c_acid, c_fe3, u_rel, feoh3_pct, ore_frac: np.ndarray,
                            porosity: np.ndarray, reactions=(True, True, True, True)) -> np.ndarray:
        """Скорость перевода урана из руды в раствор, кг/м3 в сут."""
        r = self.rate(c_acid, c_fe3, u_rel, feoh3_pct, reactions)
        return ORE_DENSITY * np.maximum(1.0 - porosity, 0.0) * ore_frac * r

    def acid_rate(self, uptake, reactions=(True, True, True, True)) -> np.ndarray:
        """Изменение концентрации кислоты в растворе, г/л в сут."""
        return -ACID_PER_U * uptake

    def uranium_rate(self, uptake) -> np.ndarray:
        """Изменение концентрации урана в растворе, мг/л в сут."""
        return 1000.0 * uptake

    def fe3_rate(self, uptake, reactions=(True, True, True, True)) -> np.ndarray:
        return -FE3_PER_U * uptake


@dataclass
class MassTransferResult:
    times: list[float] = field(default_factory=list)
    dates: list[str] = field(default_factory=list)
    acid: list[np.ndarray] = field(default_factory=list)        # г/л
    uranium: list[np.ndarray] = field(default_factory=list)     # мг/л
    residual: list[np.ndarray] = field(default_factory=list)    # кг/т
    ur_rel: list[np.ndarray] = field(default_factory=list)      # доля остаточного U
    pr_uranium: list[float] = field(default_factory=list)       # мг/л в откачиваемом растворе
    vr_acid: list[float] = field(default_factory=list)          # г/л в выщелачивающем растворе
    cumulative_u: list[float] = field(default_factory=list)     # кг извлечённого U
    cumulative_acid: list[float] = field(default_factory=list)  # кг израсходованной H2SO4
    prod_volume: list[float] = field(default_factory=list)      # м3 откачано
    inj_volume: list[float] = field(default_factory=list)       # м3 закачано
    history: list[dict] = field(default_factory=list)


def _injection_cells(wells, grid: GridSpec) -> dict[tuple[int, int, int], dict]:
    """Ячейки закачки с удельным расходом (м3/сут) на ячейку."""
    out: dict[tuple[int, int, int], dict] = {}
    for w in wells:
        if w.well_type != WellType.INJECTION:
            continue
        cells = w.screen_cells(grid)
        if not cells:
            continue
        per_cell = w.rate / len(cells)
        for c in cells:
            out[c] = {"acid": w.acid_in, "well": w.id, "rate": per_cell, "total_rate": w.rate}
    return out


def _production_cells(wells, grid: GridSpec) -> dict[tuple[int, int, int], dict]:
    """Ячейки отбора с удельным расходом (м3/сут) на ячейку."""
    out: dict[tuple[int, int, int], dict] = {}
    for w in wells:
        if w.well_type != WellType.PRODUCTION:
            continue
        cells = w.screen_cells(grid)
        if not cells:
            continue
        per_cell = abs(w.rate) / len(cells)
        for c in cells:
            out[c] = {"well": w.id, "rate": per_cell, "total_rate": abs(w.rate)}
    return out


def simulate_3d(
    grid: GridSpec,
    k: np.ndarray,
    porosity: np.ndarray,
    uranium_grade: np.ndarray,
    clay: np.ndarray,
    face_x: np.ndarray,
    face_y: np.ndarray,
    face_z: np.ndarray,
    wells,
    params: ChemicalParams,
    output_times: list[float] | None = None,
    reactions: tuple[bool, bool, bool, bool] = (True, True, True, True),
    cfl: float = 0.6,
    max_steps: int = 400_000,
    start_date: str | None = None,
    progress=None,
) -> MassTransferResult:
    """Расчёт массопереноса в 3D (адвекция + дисперсия + химические реакции).

    face_x, face_y, face_z — объёмные расходы через грани (м3/сут).
    uranium_grade — начальное содержание урана в руде, кг/т.
    """
    kin = Kinetics(params)
    nx, ny, nz = grid.shape
    shape = grid.shape

    # Начальные условия
    c_acid = np.zeros(shape)
    c_u = np.zeros(shape)
    c_fe3 = np.full(shape, params.fe3_conc)
    ore_frac = np.clip(np.asarray(uranium_grade, dtype=float) / 1000.0, 0.0, None)
    u_rel = np.ones(shape)
    u_rel0 = np.ones(shape)

    inj_cells = _injection_cells(wells, grid)
    prod_cells = _production_cells(wells, grid)

    # Начальная кислотность в растворе при наличии фактических данных
    if params.use_acid_data:
        # начальная кислотность задаётся по выщелачивающему раствору на входе
        c_acid[:] = params.acid_concentration
    else:
        c_acid[:] = params.acid_concentration

    # Начальное содержание кислоты в ячейках закачки
    for c, info in inj_cells.items():
        c_acid[c] = info["acid"]

    res = MassTransferResult()
    res.times.append(0.0)
    res.acid.append(c_acid.copy())
    res.uranium.append(c_u.copy())
    res.residual.append(uranium_grade.copy())
    res.ur_rel.append(u_rel.copy())
    res.pr_uranium.append(0.0)
    res.vr_acid.append(params.acid_concentration)
    res.cumulative_u.append(0.0)
    res.cumulative_acid.append(0.0)
    res.prod_volume.append(0.0)
    res.inj_volume.append(0.0)

    # Проводимости граней для дисперсии
    def k_face(a, b, axis):
        s = a + b
        out = np.zeros_like(s, dtype=float)
        nz = s > 0
        out[nz] = 2.0 * a[nz] * b[nz] / s[nz]
        return out

    if k.ndim == 3:
        kf_x = k_face(k[:-1, :, :], k[1:, :, :], 0)
        kf_y = k_face(k[:, :-1, :], k[:, 1:, :], 1)
        kf_z = k_face(k[:, :, :-1], k[:, :, 1:], 2)
    else:
        kf_x = kf_y = kf_z = k

    alpha_l = 0.05  # продольная дисперсия, м
    d_mol = 0.01    # молекулярная дисперсия, м2/сут
    phi = np.where(np.asarray(porosity) > 1e-6, np.asarray(porosity), 1e-6)
    cell_vol = grid.cell_volume
    pore_vol = phi * cell_vol

    # Объёмные расходы через грани (м3/сут) — из гидродинамического модуля
    Fx = np.asarray(face_x, dtype=float)
    Fy = np.asarray(face_y, dtype=float)
    Fz = np.asarray(face_z, dtype=float)

    def outflux(F, axis):
        """Суммарный исходящий расход из каждой ячейки по грани, м3/сут."""
        if axis == 0:
            out = np.zeros(shape)
            np.add.at(out, (slice(0, -1), slice(None), slice(None)), np.maximum(F, 0))
            np.add.at(out, (slice(1, None), slice(None), slice(None)), np.maximum(-F, 0))
        elif axis == 1:
            out = np.zeros(shape)
            np.add.at(out, (slice(None), slice(0, -1), slice(None)), np.maximum(F, 0))
            np.add.at(out, (slice(None), slice(1, None), slice(None)), np.maximum(-F, 0))
        else:
            out = np.zeros(shape)
            np.add.at(out, (slice(None), slice(None), slice(0, -1)), np.maximum(F, 0))
            np.add.at(out, (slice(None), slice(None), slice(1, None)), np.maximum(-F, 0))
        return out

    out_total = outflux(Fx, 0) + outflux(Fy, 1) + outflux(Fz, 2)

    # ---- ограничение устойчивости явной схемы ----
    # Учитываются как адвективная составляющая (upwind), так и
    # диффузионно-дисперсионная (центральная разностная схема).
    face_area = [(grid.dy * grid.dz), (grid.dx * grid.dz), (grid.dx * grid.dy)]
    face_dist = [grid.dx, grid.dy, grid.dz]
    diff_coef = []
    for axis, F in enumerate((Fx, Fy, Fz)):
        area = face_area[axis]
        v_face = np.abs(F) / area
        d_coef = alpha_l * v_face + d_mol
        diff_coef.append(d_coef * area / face_dist[axis])   # м3/сут на грань

    stab = out_total.copy()
    for axis, dc in enumerate(diff_coef):
        if axis == 0:
            stab[:-1, :, :] += 2.0 * dc
            stab[1:, :, :] += 2.0 * dc
        elif axis == 1:
            stab[:, :-1, :] += 2.0 * dc
            stab[:, 1:, :] += 2.0 * dc
        else:
            stab[:, :, :-1] += 2.0 * dc
            stab[:, :, 1:] += 2.0 * dc

    with np.errstate(divide="ignore", invalid="ignore"):
        dt_local = np.where(stab > 1e-12, pore_vol / stab, np.inf)
    finite = dt_local[np.isfinite(dt_local)]
    dt_stability = float(finite.min()) if finite.size else 1.0

    total_time = float(params.integration_days)
    out_times = sorted(set([float(t) for t in (output_times or [])] + [total_time]))
    out_times = [t for t in out_times if t <= total_time] or [total_time]
    if out_times[-1] < total_time:
        out_times.append(total_time)

    t = 0.0
    next_out = 0
    cum_u = 0.0
    cum_acid = 0.0
    prod_vol = 0.0
    inj_vol = 0.0
    prev_u = 0.0
    prev_vol = 0.0
    step_count = 0

    def record(time_val: float) -> None:
        nonlocal cum_u, cum_acid, prod_vol, inj_vol, prev_u, prev_vol
        # Мгновенная концентрация в откачиваемом растворе за интервал
        d_u = cum_u - prev_u
        d_v = prod_vol - prev_vol
        # Мгновенная концентрация в откачиваемом растворе за интервал.
        # M[кг]/V[м3] → г/л → ×1000 → мг/л
        pr_u = (d_u * 1000.0 / d_v) if d_v > 1e-9 else 0.0
        # Кислотность выщелачивающего раствора в рудной зоне (средняя по
        # рудным ячейкам) — показывает расход кислоты при выщелачивании
        ore_mask = ore_frac > 0
        if ore_mask.any():
            vr_a = float(np.mean(c_acid[ore_mask]))
        else:
            vr_a = float(np.mean(c_acid))
        res.times.append(float(time_val))
        res.acid.append(c_acid.copy())
        res.uranium.append(c_u.copy())
        res.residual.append(uranium_grade * u_rel)
        res.ur_rel.append(u_rel.copy())
        res.pr_uranium.append(float(pr_u))
        res.vr_acid.append(float(vr_a))
        res.cumulative_u.append(float(cum_u))
        res.cumulative_acid.append(float(cum_acid))
        res.prod_volume.append(float(prod_vol))
        res.inj_volume.append(float(inj_vol))
        prev_u = cum_u
        prev_vol = prod_vol

    while t < total_time - 1e-9 and step_count < max_steps:
        # Шаг по CFL: dt <= cfl * pore_vol / outflow
        with np.errstate(divide="ignore", invalid="ignore"):
            dt_local = np.where(stab > 1e-12, cfl * pore_vol / stab, np.inf)
        finite = dt_local[np.isfinite(dt_local)]
        dt = float(finite.min()) if finite.size else 1.0
        dt = min(dt, total_time - t, 2.0)
        if dt <= 0:
            break

        # --- адвекция и дисперсия по граням ---
        def upwind_flux(F, c_lo, c_hi, dcoef):
            """Поток растворённого вещества через грань, мг/сут.

            F — объёмный расход (м3/сут) со знаком (+ из lo в hi);
            dcoef — диффузионно-дисперсионная проводимость грани (м3/сут).
            """
            F_pos = np.maximum(F, 0.0)
            F_neg = np.minimum(F, 0.0)
            adv = F_pos * c_lo + F_neg * c_hi          # адвективная часть
            diff = -dcoef * (c_hi - c_lo)               # диффузионная часть
            return adv + diff

        d_a = np.zeros(shape)
        d_u = np.zeros(shape)
        d_fe = np.zeros(shape)

        for axis, (F, c_lo, c_hi) in enumerate([
            (Fx, c_acid[:-1, :, :], c_acid[1:, :, :]),
            (Fy, c_acid[:, :-1, :], c_acid[:, 1:, :]),
            (Fz, c_acid[:, :, :-1], c_acid[:, :, 1:]),
        ]):
            if axis == 0:
                d_a[:-1, :, :] -= upwind_flux(F, c_lo, c_hi, diff_coef[0]) * dt
                d_a[1:, :, :] += upwind_flux(F, c_lo, c_hi, diff_coef[0]) * dt
            elif axis == 1:
                d_a[:, :-1, :] -= upwind_flux(F, c_lo, c_hi, diff_coef[1]) * dt
                d_a[:, 1:, :] += upwind_flux(F, c_lo, c_hi, diff_coef[1]) * dt
            else:
                d_a[:, :, :-1] -= upwind_flux(F, c_lo, c_hi, diff_coef[2]) * dt
                d_a[:, :, 1:] += upwind_flux(F, c_lo, c_hi, diff_coef[2]) * dt

        for axis, (F, c_lo, c_hi) in enumerate([
            (Fx, c_u[:-1, :, :], c_u[1:, :, :]),
            (Fy, c_u[:, :-1, :], c_u[:, 1:, :]),
            (Fz, c_u[:, :, :-1], c_u[:, :, 1:]),
        ]):
            if axis == 0:
                d_u[:-1, :, :] -= upwind_flux(F, c_lo, c_hi, diff_coef[0]) * dt
                d_u[1:, :, :] += upwind_flux(F, c_lo, c_hi, diff_coef[0]) * dt
            elif axis == 1:
                d_u[:, :-1, :] -= upwind_flux(F, c_lo, c_hi, diff_coef[1]) * dt
                d_u[:, 1:, :] += upwind_flux(F, c_lo, c_hi, diff_coef[1]) * dt
            else:
                d_u[:, :, :-1] -= upwind_flux(F, c_lo, c_hi, diff_coef[2]) * dt
                d_u[:, :, 1:] += upwind_flux(F, c_lo, c_hi, diff_coef[2]) * dt

        for axis, (F, c_lo, c_hi) in enumerate([
            (Fx, c_fe3[:-1, :, :], c_fe3[1:, :, :]),
            (Fy, c_fe3[:, :-1, :], c_fe3[:, 1:, :]),
            (Fz, c_fe3[:, :, :-1], c_fe3[:, :, 1:]),
        ]):
            if axis == 0:
                d_fe[:-1, :, :] -= upwind_flux(F, c_lo, c_hi, diff_coef[0]) * dt
                d_fe[1:, :, :] += upwind_flux(F, c_lo, c_hi, diff_coef[0]) * dt
            elif axis == 1:
                d_fe[:, :-1, :] -= upwind_flux(F, c_lo, c_hi, diff_coef[1]) * dt
                d_fe[:, 1:, :] += upwind_flux(F, c_lo, c_hi, diff_coef[1]) * dt
            else:
                d_fe[:, :, :-1] -= upwind_flux(F, c_lo, c_hi, diff_coef[2]) * dt
                d_fe[:, :, 1:] += upwind_flux(F, c_lo, c_hi, diff_coef[2]) * dt

        # --- химические реакции ---
        uptake = kin.uranium_uptake_rate(c_acid, c_fe3, u_rel, params.feoh3_conc, ore_frac, phi, reactions)
        uptake = np.where(u_rel > 1e-8, uptake, 0.0)
        # В ячейках закачки раствор постоянно обновляется свежим реагентом —
        # реакция там не учитывается, чтобы не нарушать баланс массы урана.
        for c in inj_cells:
            uptake[c] = 0.0
        # ограничение, чтобы доля не ушла в ноль за шаг
        max_uptake = u_rel * u_rel0 * ORE_DENSITY * (1 - phi) * ore_frac / dt
        uptake = np.minimum(uptake, np.maximum(max_uptake, 0.0))

        d_a += kin.acid_rate(uptake, reactions) * dt * pore_vol
        d_u += kin.uranium_rate(uptake) * dt * pore_vol
        d_fe += kin.fe3_rate(uptake, reactions) * dt * pore_vol

        # потери урана (осаждение k3 и сорбция k4).
        # loss_mass в единицах (мг/л · м3) — как и остальные члены d_u.
        if reactions[2] or reactions[3]:
            loss_rate = np.zeros(shape)
            if reactions[2]:
                loss_rate += params.k3 * c_u
            if reactions[3]:
                loss_rate += params.k4 * c_u * np.clip(clay, 0.0, 1.0)
            loss_mass = loss_rate * dt * pore_vol
            d_u -= loss_mass
            # потерянный уран возвращается в рудную фазу (как нерастворённый):
            # loss_mass [мг/л·м3] -> кг: ×1000 (л/м3) ×1e-6 (кг/мг) = ×1e-3
            loss_kg = loss_mass * 1e-3
            ore_back = loss_kg / np.maximum(ORE_DENSITY * (1 - phi) * ore_frac * u_rel0, 1e-12)
            u_rel = np.clip(u_rel + np.where(u_rel > 0, ore_back, 0.0), 0.0, 1.0)

        # --- обновление полей ---
        c_acid = np.clip(c_acid + d_a / np.maximum(pore_vol, 1e-12), 0.0, 1e4)
        c_u = np.clip(c_u + d_u / np.maximum(pore_vol, 1e-12), 0.0, 1e6)
        c_fe3 = np.clip(c_fe3 + d_fe / np.maximum(pore_vol, 1e-12), 0.0, 1e4)

        # растворённый уран уменьшает остаточное содержание в руде
        rel_change = uptake * dt / np.maximum(ORE_DENSITY * (1 - phi) * ore_frac * u_rel0, 1e-12)
        u_rel = np.clip(u_rel - rel_change, 0.0, 1.0)

        # --- отбор продуктивного раствора ---
        step_prod_vol = 0.0
        step_u_mass = 0.0
        for c, info in prod_cells.items():
            q = info["rate"]
            step_prod_vol += q * dt
            step_u_mass += c_u[c] * 1e-3 * q * dt  # мг/л -> кг
            # извлечение вещества из ячейки
            d_mass = c_u[c] * 1e-3 * q * dt
            c_u[c] = max(0.0, c_u[c] - d_mass / max(pore_vol[c], 1e-12) * 1000.0)
            c_acid[c] = max(0.0, c_acid[c] - (c_acid[c] * 1e-3 * q * dt) / max(pore_vol[c], 1e-12) * 1000.0)
            c_fe3[c] = max(0.0, c_fe3[c] - (c_fe3[c] * 1e-3 * q * dt) / max(pore_vol[c], 1e-12) * 1000.0)
        prod_vol += step_prod_vol
        cum_u += step_u_mass
        # расход кислоты: закачка минус вынос
        for c, info in inj_cells.items():
            inj_vol += info["rate"] * dt
        cum_acid += float(np.sum(np.maximum(kin.acid_rate(uptake, reactions), 0.0) * dt * pore_vol)) / 1000.0

        # --- источники в закачивающих ячейках ---
        for c, info in inj_cells.items():
            c_acid[c] = info["acid"]
            c_u[c] = 0.0
            c_fe3[c] = params.fe3_conc

        t += dt
        step_count += 1

        if next_out < len(out_times) and t >= out_times[next_out] - 1e-9:
            record(out_times[next_out])
            next_out += 1
            if progress:
                progress(min(1.0, t / total_time), t)

    if len(res.times) <= 1 or res.times[-1] < total_time - 1e-6:
        record(total_time)

    return res


def simulate_flowline(
    points: np.ndarray,
    times: np.ndarray,
    grade: np.ndarray,
    porosity: np.ndarray,
    params: ChemicalParams,
    reactions: tuple[bool, bool, bool, bool] = (True, True, True, True),
) -> dict:
    """Расчёт концентраций вдоль линии тока (1D plug-flow по времени пребывания).

    На каждом участке линии раствор проводит время dt = times[i+1]-times[i],
    за которое происходит растворение урана по кинетическим уравнениям.
    """
    kin = Kinetics(params)
    n = len(points)
    if n < 2:
        return {"acid": [], "uranium": [], "residual": [], "times": []}

    c_acid = params.flowline_acid
    c_u = 0.0
    c_fe3 = params.fe3_conc
    u_rel = 1.0

    acid_prof = [c_acid]
    u_prof = [c_u]
    rel_prof = [u_rel]
    t_prof = [float(times[0])]
    ore_frac = np.clip(grade / 1000.0, 0.0, None)

    for i in range(n - 1):
        dt = max(0.0, float(times[i + 1]) - float(times[i]))
        if dt <= 0:
            acid_prof.append(c_acid)
            u_prof.append(c_u)
            rel_prof.append(u_rel)
            t_prof.append(float(times[i + 1]))
            continue
        # интегрирование кинетики на участке (несколько подшагов)
        steps = max(1, int(min(200, max(2, dt * 4))))
        sdt = dt / steps
        for _ in range(steps):
            ore_frac_i = ore_frac[i] if ore_frac.ndim == 1 else float(ore_frac[i])
            phi = porosity[i] if porosity.ndim == 1 else float(porosity[i])
            uptake = kin.uranium_uptake_rate(
                np.array([c_acid]), np.array([c_fe3]), np.array([u_rel]),
                np.array([params.feoh3_conc]), np.array([ore_frac_i]),
                np.array([phi]), reactions)
            up = float(uptake[0])
            # ограничение по остатку
            max_up = u_rel * ORE_DENSITY * (1 - phi) * ore_frac_i / max(sdt, 1e-9)
            up = min(up, max(0.0, max_up))
            c_acid = max(0.0, c_acid + float(kin.acid_rate(up, reactions)) * sdt)
            c_u = max(0.0, c_u + float(kin.uranium_rate(up)) * sdt)
            c_fe3 = max(0.0, c_fe3 + float(kin.fe3_rate(up, reactions)) * sdt)
            u_rel = max(0.0, u_rel - up * sdt / max(ORE_DENSITY * (1 - phi) * ore_frac_i, 1e-12))
        acid_prof.append(c_acid)
        u_prof.append(c_u)
        rel_prof.append(u_rel)
        t_prof.append(float(times[i + 1]))

    return {
        "acid": acid_prof,
        "uranium": u_prof,
        "residual": [float(r * (grade[0] if len(grade) else 0.0)) for r in rel_prof],
        "rel": rel_prof,
        "times": t_prof,
    }
