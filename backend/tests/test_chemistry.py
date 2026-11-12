"""Тесты модуля химической кинетики и массопереноса."""

from __future__ import annotations

import numpy as np
import pytest

from app import chemistry
from app.models import ChemicalParams, GridSpec, ScreenInterval, TrajectoryPoint, Well, WellType


def make_case(nx=6, ny=5, nz=3, k=2.0, phi=0.2, grade=0.8):
    grid = GridSpec(nx=nx, ny=ny, nz=nz, dx=50.0, dy=50.0, dz=10.0)
    kk = np.full(grid.shape, k)
    porosity = np.full(grid.shape, phi)
    uranium = np.full(grid.shape, grade)
    clay = np.full(grid.shape, 0.1)
    inj = Well(id="I", name="I", well_type=WellType.INJECTION,
               trajectory=[TrajectoryPoint(x=100.0, y=100.0, z=0.0)],
               screens=[ScreenInterval(top=5.0, bottom=25.0)],
               rate=500.0, acid_in=20.0)
    prod = Well(id="P", name="P", well_type=WellType.PRODUCTION,
                trajectory=[TrajectoryPoint(x=250.0, y=200.0, z=0.0)],
                screens=[ScreenInterval(top=5.0, bottom=25.0)],
                rate=-500.0)
    return grid, kk, porosity, uranium, clay, [inj, prod]


def uniform_fluxes(grid, value=100.0):
    """Простейшее поле расходов: одинаковый приток по X во все ячейки."""
    fx = np.full((grid.nx - 1, grid.ny, grid.nz), value)
    fy = np.zeros((grid.nx, grid.ny - 1, grid.nz))
    fz = np.zeros((grid.nx, grid.ny, grid.nz - 1))
    return fx, fy, fz


def test_no_reaction_keeps_concentration_uniform():
    """Без реакций концентрация остаётся однородной (перенос консервативен)."""
    grid, k, phi, grade, clay, wells = make_case()
    params = ChemicalParams(integration_days=20.0, fe3_conc=0.0, feoh3_conc=0.0,
                            beta=0.0, k1=0.0, k2=0.0, k3=0.0, k4=0.0)
    fx, fy, fz = uniform_fluxes(grid, 0.0)
    # без скважин и без реакций концентрация не меняется
    res = chemistry.simulate_3d(grid, k, phi, grade, clay, fx, fy, fz, [], params, [10.0, 20.0])
    for field in res.acid:
        assert np.allclose(field, 20.0, atol=1e-6)


def test_reaction_consumes_acid_and_releases_uranium():
    """Реакция уменьшает кислотность и увеличивает уран в растворе."""
    grid, k, phi, grade, clay, wells = make_case()
    params = ChemicalParams(integration_days=60.0, fe3_conc=1.0, feoh3_conc=0.0,
                            beta=0.0, k1=1e-3, k2=1e-4, k3=0.0, k4=0.0)
    fx, fy, fz = uniform_fluxes(grid, 0.0)
    res = chemistry.simulate_3d(grid, k, phi, grade, clay, fx, fy, fz, [], params, [60.0])
    acid_end = res.acid[-1]
    ur_end = res.uranium[-1]
    assert acid_end.mean() < 20.0            # кислота израсходована
    assert ur_end.max() > 0.0               # уран перешёл в раствор
    assert float(res.residual[-1].mean()) < float(grade.mean())  # остаток снизился


def test_temperature_increases_reaction_rate():
    """Повышение температуры ускоряет реакцию (правило Аррениуса)."""
    params = ChemicalParams(temperature=20.0, activation_energy=20.0)
    cold = chemistry.Kinetics(params)
    params_t = ChemicalParams(temperature=40.0, activation_energy=20.0)
    hot = chemistry.Kinetics(params_t)
    assert hot.temperature_factor() > cold.temperature_factor()


def test_colmatation_slows_reaction():
    """Кольматация (Fe(OH)3) замедляет выщелачивание."""
    base = ChemicalParams(k1=1e-3, feoh3_conc=0.0, beta=0.0)
    col = ChemicalParams(k1=1e-3, feoh3_conc=2.0, beta=1.0)
    kb = chemistry.Kinetics(base)
    kc = chemistry.Kinetics(col)
    args = (np.array([20.0]), np.array([1.0]), np.array([1.0]))
    assert kc.rate(*args, np.array([2.0])) < kb.rate(*args, np.array([0.0]))


def test_uranium_mass_balance():
    """Извлечённый уран не превышает начальный запас в руде."""
    grid, k, phi, grade, clay, wells = make_case()
    params = ChemicalParams(integration_days=120.0, fe3_conc=2.0, feoh3_conc=0.0,
                            beta=0.0, k1=2e-3, k2=1e-4, k3=0.0, k4=0.0)
    fx, fy, fz = uniform_fluxes(grid, 20.0)
    res = chemistry.simulate_3d(grid, k, phi, grade, clay, fx, fy, fz, wells, params, [120.0])
    pore = phi * grid.cell_volume
    ore_u = float((chemistry.ORE_DENSITY * (1 - phi) * grade / 1000.0 * grid.cell_volume).sum())
    fluid_u = float((res.uranium[-1] * pore).sum()) * 1e-6
    produced = res.cumulative_u[-1]
    assert produced <= ore_u * 0.5          # не более половины запаса за год
    assert produced + fluid_u <= ore_u * 0.6
    assert res.cumulative_u[-1] >= res.cumulative_u[0]


def test_flowline_chemistry_monotonic():
    """Вдоль линии тока кислотность убывает, уран накапливается."""
    grid, k, phi, grade, clay, wells = make_case()
    params = ChemicalParams(integration_days=100.0, fe3_conc=1.0, feoh3_conc=0.0,
                            beta=0.0, k1=5e-3, k2=1e-4, k3=0.0, k4=0.0)
    points = np.array([[100.0, 100.0, 10.0], [150.0, 150.0, 10.0],
                        [200.0, 200.0, 10.0], [250.0, 200.0, 10.0]])
    times = np.array([0.0, 10.0, 20.0, 30.0])
    prof = chemistry.simulate_flowline(points, times, np.full(4, 0.8),
                                       np.full(4, 0.2), params)
    assert prof["acid"][-1] < prof["acid"][0]
    assert prof["uranium"][-1] > prof["uranium"][0]
    assert prof["rel"][-1] < 1.0
