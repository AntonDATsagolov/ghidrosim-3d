"""Тесты гидродинамического модуля."""

from __future__ import annotations

import numpy as np
import pytest

from app.hydro import build_sources, compute_fluxes, solve
from app.models import (
    BoundaryCondition,
    GridSpec,
    ScreenInterval,
    TrajectoryPoint,
    Well,
    WellType,
)


def make_grid(nx=8, ny=6, nz=3, dx=50.0, dy=50.0, dz=10.0) -> GridSpec:
    return GridSpec(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz)


def make_well(wid, x, y, wtype, rate, z=5.0, screen=(10.0, 25.0)):
    return Well(
        id=wid, name=wid, well_type=wtype,
        trajectory=[TrajectoryPoint(x=x, y=y, z=z)],
        screens=[ScreenInterval(top=screen[0], bottom=screen[1])],
        rate=rate,
    )


def test_uniform_flow_pressure_gradient():
    """Однородный пласт: напор линейно падает от закачки к отбору."""
    grid = make_grid()
    k = np.full(grid.shape, 2.0)
    wells = [
        make_well("I1", 50.0, 50.0, WellType.INJECTION, 100.0),
        make_well("P1", 350.0, 250.0, WellType.PRODUCTION, -100.0),
    ]
    sources, _ = build_sources(grid, wells)
    res = solve(grid, k, sources, np.full(grid.shape, 0.2))
    # у закачки напор выше, у отбора — ниже
    i_in = np.unravel_index(np.argmax(sources), sources.shape)
    i_out = np.unravel_index(np.argmin(sources), sources.shape)
    assert res.head[i_in] > res.head[i_out]
    # баланс расходов соблюдён
    assert abs(res.total_in - res.total_out) < 1e-6


def test_discrete_continuity():
    """Сумма расходов по граням плюс источник равна нулю в каждой ячейке."""
    grid = make_grid()
    k = np.full(grid.shape, 1.5)
    k[2:5, 2:4, :] = 0.4
    wells = [
        make_well("I1", 75.0, 75.0, WellType.INJECTION, 300.0),
        make_well("P1", 325.0, 225.0, WellType.PRODUCTION, -300.0),
    ]
    sources, _ = build_sources(grid, wells)
    res = solve(grid, k, sources, np.full(grid.shape, 0.18), kv_ratio=0.3)
    qx, qy, qz, fx, fy, fz = compute_fluxes(grid, k, res.head, 0.3)
    inflow = np.zeros(grid.shape)
    for i in range(grid.nx):
        for j in range(grid.ny):
            for kk in range(grid.nz):
                v = 0.0
                if i > 0:
                    v += fx[i - 1, j, kk]
                if i < grid.nx - 1:
                    v -= fx[i, j, kk]
                if j > 0:
                    v += fy[i, j - 1, kk]
                if j < grid.ny - 1:
                    v -= fy[i, j, kk]
                if kk > 0:
                    v += fz[i, j, kk - 1]
                if kk < grid.nz - 1:
                    v -= fz[i, j, kk]
                inflow[i, j, kk] = v
    # приток извне равен отбору источника (знак источника обратен)
    assert np.abs(inflow + sources).max() < 1e-4 * max(1.0, abs(sources).max())


def test_dirichlet_boundary_produces_outflow():
    """Дирихле на грани: при заданном напоре из блока выходит вода."""
    grid = make_grid()
    k = np.full(grid.shape, 2.0)
    sources = np.zeros(grid.shape)
    # два грани с разным напором задают градиент и течение через блок
    faces = {"x_min": BoundaryCondition.DIRICHLET, "x_max": BoundaryCondition.DIRICHLET}
    heads = {"x_min": 10.0, "x_max": 0.0}
    res = solve(grid, k, sources, np.full(grid.shape, 0.2), faces, heads)
    assert res.boundary_out > 0
    assert res.boundary_in > 0
    assert abs(res.residual) < 1e-6


def test_neumann_blocks_flow():
    """При непроницаемых границах и отсутствии скважин потока нет."""
    grid = make_grid()
    k = np.full(grid.shape, 2.0)
    res = solve(grid, k, np.zeros(grid.shape), np.full(grid.shape, 0.2))
    assert np.allclose(res.speed, 0.0, atol=1e-9)
    assert abs(res.head.max() - res.head.min()) < 1e-6


def test_anisotropy_reduces_vertical_flow():
    """При k_z << k_xy вертикальные расходы меньше горизонтальных."""
    grid = make_grid()
    k = np.full(grid.shape, 2.0)
    wells = [make_well("I1", 100.0, 100.0, WellType.INJECTION, 200.0)]
    sources, _ = build_sources(grid, wells)
    res = solve(grid, k, sources, np.full(grid.shape, 0.2), kv_ratio=0.05)
    assert np.abs(res.darcy_z).max() < np.abs(res.darcy_x).max()


def test_screen_cells_cover_interval():
    """Ячейки фильтра выбираются по всему интервалу перфорации."""
    grid = make_grid(nz=5, dz=10.0)
    w = make_well("W", 125.0, 125.0, WellType.INJECTION, 100.0, screen=(15.0, 45.0))
    cells = w.screen_cells(grid)
    # центры слоёв: 5, 15, 25, 35, 45 (шаг 10) -> в интервал попадают 15,25,35,45
    assert len(cells) == 4


def test_source_distribution_matches_well_rate():
    grid = make_grid()
    wells = [make_well("I1", 125.0, 125.0, WellType.INJECTION, 400.0),
             make_well("P1", 325.0, 225.0, WellType.PRODUCTION, -400.0)]
    sources, cells_map = build_sources(grid, wells)
    assert abs(sources.sum()) < 1e-9
    assert sources.max() > 0 and sources.min() < 0
