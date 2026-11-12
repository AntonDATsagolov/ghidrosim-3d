"""Трассировка линий тока (стримлайнов) по полю скоростей.

Используется для:
  * визуализации направления движения раствора между скважинами;
  * расчёта времени пребывания раствора по линии тока;
  * расчёта массопереноса вдоль линий тока;
  * выявления застойных зон (где скорость ниже порога).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .models import GridSpec


@dataclass
class Streamline:
    points: np.ndarray          # (n, 3) — координаты X, Y, Z, м
    speeds: np.ndarray          # (n,) — скорость, м/сут
    times: np.ndarray           # (n,) — время от начала, сут
    well_in: str = ""
    well_out: str = ""
    stagnant: bool = False
    residence: float = 0.0      # полное время прохождения линии, сут

    def to_dict(self) -> dict:
        return {
            "points": self.points.tolist(),
            "speeds": self.speeds.tolist(),
            "times": self.times.tolist(),
            "well_in": self.well_in,
            "well_out": self.well_out,
            "stagnant": self.stagnant,
            "residence": self.residence,
        }


def _trilinear(field: np.ndarray, grid: GridSpec, p: np.ndarray) -> np.ndarray:
    """Трилинейная интерполяция поля в точках p (N,3) -> (N,)."""
    nx, ny, nz = grid.shape
    fx = (p[:, 0] - grid.x0) / grid.dx - 0.5
    fy = (p[:, 1] - grid.y0) / grid.dy - 0.5
    fz = (p[:, 2] - grid.z0) / grid.dz - 0.5
    fx = np.clip(fx, 0, nx - 1.000001)
    fy = np.clip(fy, 0, ny - 1.000001)
    fz = np.clip(fz, 0, nz - 1.000001)
    i0 = np.floor(fx).astype(int); j0 = np.floor(fy).astype(int); k0 = np.floor(fz).astype(int)
    i1 = np.minimum(i0 + 1, nx - 1); j1 = np.minimum(j0 + 1, ny - 1); k1 = np.minimum(k0 + 1, nz - 1)
    tx = fx - i0; ty = fy - j0; tz = fz - k0

    def g(ii, jj, kk):
        return field[ii, jj, kk]

    c00 = g(i0, j0, k0) * (1 - tx) + g(i1, j0, k0) * tx
    c10 = g(i0, j1, k0) * (1 - tx) + g(i1, j1, k0) * tx
    c01 = g(i0, j0, k1) * (1 - tx) + g(i1, j0, k1) * tx
    c11 = g(i0, j1, k1) * (1 - tx) + g(i1, j1, k1) * tx
    c0 = c00 * (1 - ty) + c10 * ty
    c1 = c01 * (1 - ty) + c11 * ty
    return c0 * (1 - tz) + c1 * tz


def _velocity(vx, vy, vz, grid: GridSpec, p: np.ndarray) -> np.ndarray:
    return np.stack(
        [
            _trilinear(vx, grid, p),
            _trilinear(vy, grid, p),
            _trilinear(vz, grid, p),
        ],
        axis=1,
    )


def trace(
    grid: GridSpec,
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    seeds: np.ndarray,
    seed_labels: list[str] | None = None,
    stop_cells: dict[tuple[int, int, int], str] | None = None,
    max_steps: int = 4000,
    stagnation_threshold: float = 0.05,
    max_length: float | None = None,
    max_residence: float = 3650.0,
) -> list[Streamline]:
    """Трассировка линий тока методом РК4 по полю скоростей.

    seeds — массив (N, 3) начальных точек.
    stop_cells — словарь {индекс ячейки: идентификатор скважины-приёма}.
    """
    seeds = np.atleast_2d(np.asarray(seeds, dtype=float))
    if len(seeds) == 0:
        return []
    seed_labels = seed_labels or [""] * len(seeds)
    stop_cells = stop_cells or {}

    n = len(seeds)
    step = 0.35 * min(grid.dx, grid.dy, grid.dz)
    if max_length:
        max_steps = int(min(max_steps, max_length / max(step, 1e-9)))

    P = seeds.copy()
    alive = np.ones(n, dtype=bool)
    out_well = np.array([""] * n, dtype=object)
    stagnant = np.zeros(n, dtype=bool)

    paths: list[list[np.ndarray]] = [[p] for p in P]
    speeds_l: list[list[float]] = [[0.0] for _ in range(n)]
    times_l: list[list[float]] = [[0.0] for _ in range(n)]
    cur_t = np.zeros(n)
    cur_s = np.zeros(n)

    # начальные скорости
    v = _velocity(vx, vy, vz, grid, P)
    sp0 = np.sqrt((v**2).sum(axis=1))
    for i in range(n):
        speeds_l[i][0] = float(sp0[i])

    for it in range(max_steps):
        if not alive.any():
            break
        idx_alive = np.where(alive)[0]
        Pa = P[idx_alive]
        va = v[idx_alive]

        norm = np.sqrt((va**2).sum(axis=1))
        norm = np.where(norm < 1e-14, 1e-14, norm)
        dirv = va / norm[:, None]

        k1 = dirv
        k2 = dirv + 0.0
        v_mid = _velocity(vx, vy, vz, grid, Pa + 0.5 * step * k1)
        n2 = np.sqrt((v_mid**2).sum(axis=1))
        n2 = np.where(n2 < 1e-14, 1e-14, n2)
        k2 = v_mid / n2[:, None]
        v_mid2 = _velocity(vx, vy, vz, grid, Pa + 0.5 * step * k2)
        n3 = np.sqrt((v_mid2**2).sum(axis=1))
        n3 = np.where(n3 < 1e-14, 1e-14, n3)
        k3 = v_mid2 / n3[:, None]
        v_end = _velocity(vx, vy, vz, grid, Pa + step * k3)
        n4 = np.sqrt((v_end**2).sum(axis=1))
        n4 = np.where(n4 < 1e-14, 1e-14, n4)
        k4 = v_end / n4[:, None]

        newP = Pa + (step / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        vnew = _velocity(vx, vy, vz, grid, newP)
        snew = np.sqrt((vnew**2).sum(axis=1))

        # обновляем состояние
        P[idx_alive] = newP
        v[idx_alive] = vnew
        cur_s[idx_alive] = snew
        dt = step / np.maximum(norm, 1e-9)
        cur_t[idx_alive] += dt

        for local, gi in enumerate(idx_alive):
            paths[gi].append(newP[local].copy())
            speeds_l[gi].append(float(snew[local]))
            times_l[gi].append(float(cur_t[gi]))

        # остановки
        for local, gi in enumerate(idx_alive):
            p = newP[local]
            # выход за пределы блока
            if (p[0] < grid.x0 - 0.5 * grid.dx or p[0] > grid.x0 + grid.nx * grid.dx + 0.5 * grid.dx or
                    p[1] < grid.y0 - 0.5 * grid.dy or p[1] > grid.y0 + grid.ny * grid.dy + 0.5 * grid.dy or
                    p[2] < grid.z0 - 0.5 * grid.dz or p[2] > grid.z0 + grid.nz * grid.dz + 0.5 * grid.dz):
                alive[gi] = False
                continue
            if snew[local] < stagnation_threshold:
                stagnant[gi] = True
                alive[gi] = False
                continue
            if cur_t[gi] > max_residence:
                alive[gi] = False
                continue
            ci = int(np.clip(round((p[0] - grid.x0) / grid.dx - 0.5), 0, grid.nx - 1))
            cj = int(np.clip(round((p[1] - grid.y0) / grid.dy - 0.5), 0, grid.ny - 1))
            ck = int(np.clip(round((p[2] - grid.z0) / grid.dz - 0.5), 0, grid.nz - 1))
            if (ci, cj, ck) in stop_cells:
                out_well[gi] = stop_cells[(ci, cj, ck)]
                alive[gi] = False

    result: list[Streamline] = []
    for i in range(n):
        pts = np.array(paths[i])
        if len(pts) < 2:
            continue
        result.append(
            Streamline(
                points=pts,
                speeds=np.array(speeds_l[i]),
                times=np.array(times_l[i]),
                well_in=seed_labels[i] if i < len(seed_labels) else "",
                well_out=str(out_well[i]),
                stagnant=bool(stagnant[i]),
                residence=float(cur_t[i]),
            )
        )
    return result


def seeds_from_wells(wells, grid: GridSpec, well_types=("injection",)) -> tuple[np.ndarray, list[str]]:
    """Начальные точки линий тока — центры ячеек фильтров закачивающих скважин."""
    pts: list[np.ndarray] = []
    labels: list[str] = []
    for w in wells:
        if w.well_type.value not in well_types:
            continue
        for (i, j, k) in w.screen_cells(grid):
            x = grid.x0 + (i + 0.5) * grid.dx
            y = grid.y0 + (j + 0.5) * grid.dy
            z = grid.z0 + (k + 0.5) * grid.dz
            pts.append([x, y, z])
            labels.append(w.id)
    if not pts:
        return np.empty((0, 3)), []
    return np.array(pts), labels


def stop_cells_from_wells(wells, grid: GridSpec, well_types=("production",)) -> dict[tuple[int, int, int], str]:
    """Ячейки фильтров приёмных скважин — точки остановки линий тока."""
    stops: dict[tuple[int, int, int], str] = {}
    for w in wells:
        if w.well_type.value not in well_types:
            continue
        for c in w.screen_cells(grid):
            stops[c] = w.id
    return stops


def stagnation_map(grid: GridSpec, speed: np.ndarray, threshold: float) -> np.ndarray:
    """Карта застойных зон (доля ячеек с малой скоростью в прорезях)."""
    stag = speed < threshold
    out = np.zeros(grid.shape)
    out[stag] = 1.0
    return out
