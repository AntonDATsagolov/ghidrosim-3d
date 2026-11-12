"""Вспомогательные геометрические операции: интерполяция полей в точку,
построение интервалов фильтрации, работа с траекторией скважин."""

from __future__ import annotations

import numpy as np

from .models import GridSpec, ScreenInterval, TrajectoryPoint, Well


def make_screen(d: dict) -> ScreenInterval:
    return ScreenInterval(top=float(d["top"]), bottom=float(d["bottom"]))


def sample_field(field: np.ndarray, grid: GridSpec, p) -> float:
    """Значение поля в точке (трилинейная интерполяция, с отсечкой по границам)."""
    nx, ny, nz = grid.shape
    fx = (float(p[0]) - grid.x0) / grid.dx - 0.5
    fy = (float(p[1]) - grid.y0) / grid.dy - 0.5
    fz = (float(p[2]) - grid.z0) / grid.dz - 0.5
    fx = np.clip(fx, 0, nx - 1.000001)
    fy = np.clip(fy, 0, ny - 1.000001)
    fz = np.clip(fz, 0, nz - 1.000001)
    i0, j0, k0 = int(fx), int(fy), int(fz)
    i1, j1, k1 = min(i0 + 1, nx - 1), min(j0 + 1, ny - 1), min(k0 + 1, nz - 1)
    tx, ty, tz = fx - i0, fy - j0, fz - k0

    def g(i, j, k):
        return float(field[i, j, k])

    c00 = g(i0, j0, k0) * (1 - tx) + g(i1, j0, k0) * tx
    c10 = g(i0, j1, k0) * (1 - tx) + g(i1, j1, k0) * tx
    c01 = g(i0, j0, k1) * (1 - tx) + g(i1, j0, k1) * tx
    c11 = g(i0, j1, k1) * (1 - tx) + g(i1, j1, k1) * tx
    c0 = c00 * (1 - ty) + c10 * ty
    c1 = c01 * (1 - ty) + c11 * ty
    return c0 * (1 - tz) + c1 * tz


def sample_fields(fields: dict[str, np.ndarray], grid: GridSpec, p) -> dict[str, float]:
    return {k: sample_field(v, grid, p) for k, v in fields.items()}


def well_path(well: Well) -> list[TrajectoryPoint]:
    if well.trajectory:
        return list(well.trajectory)
    return []


def trajectory_length(points: list[TrajectoryPoint]) -> float:
    if len(points) < 2:
        return 0.0
    total = 0.0
    for a, b in zip(points[:-1], points[1:]):
        total += float(np.sqrt((b.x - a.x) ** 2 + (b.y - a.y) ** 2 + (b.z - a.z) ** 2))
    return total


def auto_grid_bounds(well_points: list[tuple[float, float, float]],
                     margin: float = 50.0,
                     cell: float = 25.0,
                     max_cells: int = 400_000) -> GridSpec:
    """Подбор сетки по облаку точек скважин."""
    if not well_points:
        return GridSpec()
    arr = np.array(well_points, dtype=float)
    x0, x1 = arr[:, 0].min() - margin, arr[:, 0].max() + margin
    y0, y1 = arr[:, 1].min() - margin, arr[:, 1].max() + margin
    z0, z1 = arr[:, 2].min() - margin, arr[:, 2].max() + margin
    nx = max(4, int(np.ceil((x1 - x0) / cell)))
    ny = max(4, int(np.ceil((y1 - y0) / cell)))
    nz = max(2, int(np.ceil((z1 - z0) / cell)))
    while nx * ny * nz > max_cells:
        cell *= 1.3
        nx = max(4, int(np.ceil((xx := (x1 - x0)) / cell)))
        ny = max(4, int(np.ceil((y1 - y0) / cell)))
        nz = max(2, int(np.ceil((z1 - z0) / cell)))
    return GridSpec(x0=x0, y0=y0, z0=z0, nx=nx, ny=ny, nz=nz,
                    dx=cell, dy=cell, dz=cell)
