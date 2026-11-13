"""Геостатистика: обратные взвешенные расстояния (IDW) и ординарный кригинг.

Используется для восстановления свойств пласта (коэффициент фильтрации,
пористость, содержание урана, литология) в узлах расчётной сетки по
скважинным данным. Реализованы модели вариограммы: сферическая,
экспоненциальная, гауссова и линейная.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.spatial import cKDTree


@dataclass
class Variogram:
    """Модель вариограммы: γ(h) = nugget + sill * f(h/range)."""

    model: str = "spherical"       # spherical | exponential | gaussian | linear
    sill: float = 1.0
    range_: float = 100.0
    nugget: float = 0.0

    def evaluate(self, h: np.ndarray | float) -> np.ndarray:
        h = np.asarray(h, dtype=float)
        r = max(self.range_, 1e-9)
        if self.model == "spherical":
            hr = np.clip(h / r, 0.0, 1.0)
            f = 1.5 * hr - 0.5 * hr**3
        elif self.model == "exponential":
            f = 1.0 - np.exp(-3.0 * h / r)
        elif self.model == "gaussian":
            f = 1.0 - np.exp(-3.0 * (h / r) ** 2)
        elif self.model == "linear":
            # линейная модель с насыщением на силлете
            f = np.clip(h / r, 0.0, 1.0)
        else:
            raise ValueError(f"Неизвестная модель вариограммы: {self.model}")
        return self.nugget + self.sill * f


@dataclass
class GeoSettings:
    """Настройки интерполяции."""

    method: str = "kriging"   # kriging | idw
    variogram: Variogram = None  # type: ignore[assignment]
    idw_power: float = 2.0
    idw_radius: float = 0.0     # 0 — без ограничения радиуса
    anisotropy: tuple[float, float, float] = (1.0, 1.0, 1.0)
    max_neighbors: int = 12     # для IDW
    nugget: float = 0.0

    def __post_init__(self) -> None:
        if self.variogram is None:
            self.variogram = Variogram()


def _apply_anisotropy(d: np.ndarray, aniso: np.ndarray) -> np.ndarray:
    """Приведение расстояния к нормированному виду с учётом анизотропии."""
    return np.sqrt(np.sum((d / aniso) ** 2, axis=-1))


def idw(
    points: np.ndarray,
    values: np.ndarray,
    targets: np.ndarray,
    settings: GeoSettings,
) -> np.ndarray:
    """Интерполяция методом обратных взвешенных расстояний."""
    points = np.asarray(points, dtype=float)
    values = np.asarray(values, dtype=float)
    targets = np.atleast_2d(np.asarray(targets, dtype=float))
    if len(points) == 0:
        return np.full(targets.shape[0], np.nan)

    aniso = np.asarray(settings.anisotropy, dtype=float)
    aniso = np.where(aniso > 0, aniso, 1.0)

    tree = cKDTree(points / aniso)
    radius = settings.idw_radius
    if radius > 0:
        radius = radius / float(np.min(aniso))
        neigh = tree.query_ball_point(targets / aniso, r=radius, workers=-1)
    else:
        k = min(settings.max_neighbors, len(points))
        _, neigh = tree.query(targets / aniso, k=k, workers=-1)
        neigh = np.atleast_2d(neigh)

    out = np.empty(targets.shape[0], dtype=float)
    p = settings.idw_power
    for i, idx in enumerate(neigh):
        idx = np.atleast_1d(np.asarray(idx, dtype=int))
        idx = idx[idx < len(points)]
        if idx.size == 0:
            # нет соседей — ищем ближайший в глобальном масштабе
            d_all = _apply_anisotropy(targets[i] - points, aniso)
            j = int(np.argmin(d_all))
            out[i] = values[j]
            continue
        d = _apply_anisotropy(targets[i] - points[idx], aniso)
        w = 1.0 / np.maximum(d, 1e-6) ** p
        out[i] = float(np.sum(w * values[idx]) / np.sum(w))
    return out


def ordinary_kriging(
    points: np.ndarray,
    values: np.ndarray,
    targets: np.ndarray,
    settings: GeoSettings,
) -> np.ndarray:
    """Ординарный кригинг с заданной вариограммой.

    Система: [K  1][w]   [z]
             [1ᵀ 0][µ] = [1]
    Решается методом Холецкого с малым регуляризующим добавлением к диагонали.
    """
    points = np.asarray(points, dtype=float)
    values = np.asarray(values, dtype=float)
    targets = np.atleast_2d(np.asarray(targets, dtype=float))
    n = len(points)
    if n == 0:
        return np.full(targets.shape[0], np.nan)
    if n == 1:
        return np.full(targets.shape[0], float(values[0]))

    aniso = np.asarray(settings.anisotropy, dtype=float)
    aniso = np.where(aniso > 0, aniso, 1.0)
    vg = settings.variogram

    # Матрица вариограмм между исходными точками — считается один раз.
    # Матрица вариограмм условно положительно определена; ординарный кригинг
    # решается расширенной системой, поэтому спектральной проекции здесь нет.
    dist_pp = _apply_anisotropy(points[:, None, :] - points[None, :, :], aniso)
    K_full = vg.evaluate(dist_pp)
    tree = cKDTree(points / aniso)
    # Для кригинга используем локальную выборку ближайших соседей.
    # Порядок выбора фиксирован (устойчивая сортировка по расстоянию) —
    # так результат воспроизводим и совпадает с портом на JavaScript.
    max_pts = 16
    out = np.empty(targets.shape[0], dtype=float)
    if n <= max_pts:
        neigh_all = np.tile(np.arange(n), (targets.shape[0], 1))
    else:
        diff = (targets[:, None, :] - points[None, :, :]) / aniso
        dist = np.sqrt((diff ** 2).sum(axis=2))
        neigh_all = np.argsort(dist, axis=1, kind="stable")[:, :max_pts]

    for i in range(targets.shape[0]):
        idx = np.atleast_1d(neigh_all[i]).astype(int)
        Ksub = K_full[np.ix_(idx, idx)]
        vals = values[idx]
        # Регуляризация диагонали (нуджет) для устойчивости Холецкого
        nug = max(vg.nugget, 1e-8) + max(abs(np.trace(Ksub)) / max(len(idx), 1), 1e-8) * 1e-6
        Ksub = Ksub.copy()
        np.fill_diagonal(Ksub, Ksub.diagonal() + nug)

        Kaug = np.zeros((len(idx) + 1, len(idx) + 1), dtype=float)
        Kaug[: len(idx), : len(idx)] = Ksub
        Kaug[len(idx), : len(idx)] = 1.0
        Kaug[: len(idx), len(idx)] = 1.0

        # Система ординарного кригинга в форме вариограммы:
        #   [ Γ   1 ] [λ]   [γ(h_i0)]
        #   [ 1ᵀ  0 ] [μ] = [  1   ]
        dist = _apply_anisotropy(targets[i] - points[idx], aniso)
        kvec = vg.evaluate(dist)
        b = np.concatenate([kvec, [1.0]])
        try:
            lu, piv = lu_factor(Kaug)
            sol = lu_solve((lu, piv), b)
            if not np.all(np.isfinite(sol)):
                raise ValueError('non-finite')
        except Exception:
            sol = np.linalg.lstsq(Kaug, b, rcond=None)[0]
        w = sol[: len(vals)]
        out[i] = float(np.sum(w * vals))
    return out


def _conditional_positive_semidefinite(K: np.ndarray) -> np.ndarray:
    """Проекция матрицы на класс условно положительно определённых (вариограмма)."""
    try:
        w, Q = np.linalg.eigh(K)
        if w.min() < -1e-8 * max(abs(w).max(), 1.0):
            w = np.clip(w, 0.0, None)
            K = Q @ np.diag(w) @ Q.T
            K = 0.5 * (K + K.T)
    except np.linalg.LinAlgError:
        pass
    return K


def interpolate(
    points: np.ndarray,
    values: np.ndarray,
    targets: np.ndarray,
    settings: GeoSettings,
) -> np.ndarray:
    """Точка входа: выбор метода интерполяции."""
    points = np.asarray(points, dtype=float)
    values = np.asarray(values, dtype=float)
    if points.size == 0:
        return np.full(np.atleast_2d(targets).shape[0], np.nan)
    # Отбрасываем точки с пропусками
    mask = np.isfinite(values)
    points, values = points[mask], values[mask]
    if len(points) == 0:
        return np.full(np.atleast_2d(targets).shape[0], np.nan)
    if settings.method == "idw":
        return idw(points, values, targets, settings)
    return ordinary_kriging(points, values, targets, settings)


def nearest_lithology(
    points: np.ndarray,
    labels: list[str],
    targets: np.ndarray,
) -> np.ndarray:
    """Ближайший сосед для категориальных данных (литология)."""
    points = np.asarray(points, dtype=float)
    targets = np.atleast_2d(np.asarray(targets, dtype=float))
    if len(points) == 0:
        return np.full(targets.shape[0], "", dtype=object)
    tree = cKDTree(points)
    _, idx = tree.query(targets, k=1, workers=-1)
    return np.array([labels[int(i)] for i in np.atleast_1d(idx)], dtype=object)
