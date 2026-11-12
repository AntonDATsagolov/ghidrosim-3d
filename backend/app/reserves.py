"""Подсчёт запасов урана в закисляемом объёме.

В соответствии с методикой исходного ПО выполняется расчёт:
  * суммарный запас урана в ячейках, где кислотность выше заданного порога;
  * площадь блока с учётом пользовательского контура;
  * средневзвешенное содержание урана;
  * средняя мощность рудного пласта (с учётом баланса);
  * средневзвешенный метропроцент.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .models import GridSpec


@dataclass
class ReserveResult:
    reserve_kg: float = 0.0              # суммарный запас урана, кг
    reserve_t: float = 0.0               # то же в тоннах
    acidified_volume: float = 0.0        # объём закисляемых ячеек, м3
    acidified_cells: int = 0
    area_m2: float = 0.0                 # площадь в плане (м2)
    area_ha: float = 0.0
    avg_grade: float = 0.0               # средневзвешенное содержание, кг/т
    avg_grade_pct: float = 0.0
    avg_thickness: float = 0.0           # средняя мощность, м
    weighted_metre_percent: float = 0.0  # средневзвешенный метропроцент
    min_grade: float = 0.0
    max_grade: float = 0.0
    total_ore_t: float = 0.0             # общий запас руды в контуре, т
    contour_area_m2: float = 0.0
    contour_volume: float = 0.0
    by_layer: list[dict] = field(default_factory=list)


def calculate(
    grid: GridSpec,
    grade: np.ndarray,
    acid: np.ndarray,
    acid_threshold: float = 5.0,
    contour_mask: np.ndarray | None = None,
    porosity: np.ndarray | None = None,
    density: float = 2200.0,
) -> ReserveResult:
    """Подсчёт запасов в закисляемом объёме (и в пользовательском контуре).

    grade — содержание урана в руде, кг/т (текущее, с учётом отработки).
    acid — концентрация кислоты в растворе, г/л.
    contour_mask — булева маска (nx,ny,nz) ячеек внутри заданного контура.
    """
    nx, ny, nz = grid.shape
    grade = np.asarray(grade, dtype=float)
    acid = np.asarray(acid, dtype=float)
    res = ReserveResult()

    if contour_mask is not None:
        contour_mask = np.asarray(contour_mask, dtype=bool)
        if contour_mask.shape != grid.shape:
            contour_mask = np.ones(grid.shape, dtype=bool)
    else:
        contour_mask = np.ones(grid.shape, dtype=bool)

    # Закисляемые ячейки: кислотность выше порога
    acidified = (acid >= acid_threshold) & contour_mask
    res.acidified_cells = int(np.count_nonzero(acidified))
    res.acidified_volume = res.acidified_cells * grid.cell_volume

    # Масса руды в контуре и в закисляемом объёме
    cell_ore_t = density * grid.cell_volume / 1000.0  # т руды в ячейке
    res.total_ore_t = float(np.count_nonzero(contour_mask) * cell_ore_t)

    if res.acidified_cells == 0:
        return res

    g = np.where(acidified, grade, 0.0)
    w = np.where(acidified, cell_ore_t, 0.0)
    total_w = float(w.sum())
    res.reserve_kg = float((g * w).sum())          # кг U
    res.reserve_t = res.reserve_kg / 1000.0

    # Средневзвешенное содержание, кг/т
    res.avg_grade = float((g * w).sum() / total_w) if total_w > 0 else 0.0
    res.avg_grade_pct = res.avg_grade / 10.0      # 1% = 10 кг/т
    res.min_grade = float(g[acidified].min())
    res.max_grade = float(g[acidified].max())

    # Мощность: доля ячеек руды в колонке, средняя по площади контура
    # (по оси z — пласт, по оси x,y — площадь)
    active_z = np.zeros((nx, ny))
    ore_thk = np.zeros((nx, ny))
    for i in range(nx):
        for j in range(ny):
            col = acidified[i, j, :]
            n_act = int(np.count_nonzero(col))
            if n_act:
                active_z[i, j] = n_act * grid.dz
                ore_thk[i, j] = float(grade[i, j, col].mean())
    res.avg_thickness = float(active_z.mean()) if active_z.size else 0.0

    # Средневзвешенный метропроцент: sum(мощность_i * содержание_i) / sum(мощность_i)
    mask = active_z > 0
    if np.any(mask):
        num = float((ore_thk[mask] * active_z[mask]).sum())
        den = float(active_z[mask].sum())
        res.weighted_metre_percent = num / den if den > 0 else 0.0
    else:
        res.weighted_metre_percent = 0.0

    # Площадь в плане
    cell_area = grid.dx * grid.dy
    cols_with_ore = int(np.count_nonzero(mask))
    res.area_m2 = cols_with_ore * cell_area
    res.area_ha = res.area_m2 / 10000.0
    res.contour_area_m2 = int(np.count_nonzero(contour_mask.any(axis=2))) * cell_area
    res.contour_volume = float(np.count_nonzero(contour_mask) * grid.cell_volume)

    # Разбивка по слоям (для отчёта)
    for k in range(nz):
        sl = acidified[:, :, k]
        n = int(np.count_nonzero(sl))
        if n == 0:
            continue
        gg = grade[:, :, k][sl]
        res.by_layer.append({
            "layer": k + 1,
            "z_top": grid.z0 + k * grid.dz,
            "z_bottom": grid.z0 + (k + 1) * grid.dz,
            "cells": n,
            "volume": n * grid.cell_volume,
            "reserve_kg": float((gg * cell_ore_t).sum()),
            "avg_grade": float(gg.mean()),
            "avg_acid": float(acid[:, :, k][sl].mean()),
        })

    return res


def contour_from_polygon(grid: GridSpec, polygon: list[list[float]], z_range: tuple[float, float] | None = None) -> np.ndarray:
    """Построение булевой маски контура по полигону в плане (X, Y)."""
    mask = np.zeros(grid.shape, dtype=bool)
    if not polygon or len(polygon) < 3:
        return mask
    xc, yc, zc = grid.centers()
    px = np.array([p[0] for p in polygon], dtype=float)
    py = np.array([p[1] for p in polygon], dtype=float)
    inside = _point_in_polygon(xc, yc, px, py)
    if z_range:
        inside = inside & (zc >= z_range[0]) & (zc <= z_range[1])
    return inside


def _point_in_polygon(x: np.ndarray, y: np.ndarray, px: np.ndarray, py: np.ndarray) -> np.ndarray:
    """Тест принадлежности точек полигону (алгоритм «луч»)."""
    n = len(px)
    inside = np.zeros(x.shape, dtype=bool)
    j = n - 1
    for i in range(n):
        cond = ((py[i] > y) != (py[j] > y)) & (
            x < (px[j] - px[i]) * (y - py[i]) / (py[j] - py[i] + 1e-30) + px[i]
        )
        inside ^= cond
        j = i
    return inside
