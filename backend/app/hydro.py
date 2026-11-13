"""Гидродинамический модуль: стационарная фильтрация Дарси в 3D.

Решается уравнение
        div( k * grad(h) ) + q = 0
где h — напор (м), k — коэффициент фильтрации (м/сут), q — объёмный
расход (м3/сут) с учётом знака (закачка +, отбор −).

Дискретизация — конечные объёмы по граням с использованием
проводимости между соседними ячейками (гармоническое среднее для k).

Граничные условия задаются по каждой из шести граней расчётного блока:
    * Неймана (по умолчанию) — непроницаемость (qv = 0);
    * Дирихле — заданный напор на границе (отток из блока / приток
      региональных грунтовых вод).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .models import BoundaryCondition, GridSpec

FACE_AXES = {
    "x_min": (0, 0),
    "x_max": (0, -1),
    "y_min": (1, 0),
    "y_max": (1, -1),
    "z_min": (2, 0),
    "z_max": (2, -1),
}


@dataclass
class HydroResult:
    head: np.ndarray            # напор, м (nx,ny,nz)
    vx: np.ndarray              # скорость по X, м/сут
    vy: np.ndarray
    vz: np.ndarray
    speed: np.ndarray           # модуль скорости, м/сут
    darcy_x: np.ndarray         # объёмный расход Дарси, м2/сут
    darcy_y: np.ndarray
    darcy_z: np.ndarray
    face_x: np.ndarray          # объёмный расход через грани, м3/сут
    face_y: np.ndarray
    face_z: np.ndarray
    residual: float             # невязка баланса, м3/сут
    total_in: float             # суммарная закачка, м3/сут
    total_out: float            # суммарный отбор, м3/сут
    boundary_in: float          # приток через границы, м3/сут
    boundary_out: float         # отток через границы, м3/сут
    iterations: int = 0
    solver: str = "spsolve"


def _harmonic(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    s = a + b
    out = np.zeros_like(s, dtype=float)
    nz = s > 0
    out[nz] = 2.0 * a[nz] * b[nz] / s[nz]
    return out


def _face_conductance(grid: GridSpec, k: np.ndarray, axis: int, kv_ratio: float = 1.0) -> np.ndarray:
    """Проводимость между ячейками по оси axis. Массив формы (nx,ny,nz) с нулевой
    гранью по axis, т.е. корректный размер: axis=0 -> (nx-1,ny,nz).

    kv_ratio — коэффициент анизотропии по вертикали (k_z / k_xy).
    """
    if axis == 0:
        kf = _harmonic(k[:-1, :, :], k[1:, :, :])
        return kf * (grid.dy * grid.dz) / grid.dx
    if axis == 1:
        kf = _harmonic(k[:, :-1, :], k[:, 1:, :])
        return kf * (grid.dx * grid.dz) / grid.dy
    kf = _harmonic(k[:, :, :-1], k[:, :, 1:])
    return kf * kv_ratio * (grid.dx * grid.dy) / grid.dz


def build_system(
    grid: GridSpec,
    k: np.ndarray,
    sources: np.ndarray,
    bc_faces: dict[str, BoundaryCondition] | None = None,
    bc_heads: dict[str, float] | None = None,
    kv_ratio: float = 1.0,
) -> tuple[sp.csr_matrix, np.ndarray, dict]:
    """Сборка разреженной системы для напора."""
    nx, ny, nz = grid.shape
    n = grid.n_cells
    bc_faces = bc_faces or {}
    bc_heads = bc_heads or {}

    idx = np.arange(n).reshape(nx, ny, nz)
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    data: list[np.ndarray] = []
    rhs = sources.reshape(-1).astype(float).copy()
    meta: dict = {"has_dirichlet": False}

    diag = np.zeros(n, dtype=float)
    has_dirichlet = False

    for axis in range(3):
        cond = _face_conductance(grid, k, axis, kv_ratio)  # проводимости граней
        if axis == 0:
            i0, i1 = idx[:-1, :, :], idx[1:, :, :]
            shape = (nx - 1, ny, nz)
        elif axis == 1:
            i0, i1 = idx[:, :-1, :], idx[:, 1:, :]
            shape = (nx, ny - 1, nz)
        else:
            i0, i1 = idx[:, :, :-1], idx[:, :, 1:]
            shape = (nx, ny, nz - 1)

        c = cond.ravel()
        a = i0.ravel()
        b = i1.ravel()
        rows.append(a); cols.append(b); data.append(-c)
        rows.append(b); cols.append(a); data.append(-c)
        np.add.at(diag, a, c)
        np.add.at(diag, b, c)

    # Граничные условия Дирихле на внешних гранях блока
    for face, (axis, side) in FACE_AXES.items():
        cond_bc = bc_faces.get(face, BoundaryCondition.NEUMANN)
        if cond_bc != BoundaryCondition.DIRICHLET:
            continue
        has_dirichlet = True
        h_d = float(bc_heads.get(face, 0.0))
        # Ячейки, примыкающие к грани
        if axis == 0:
            sl = [slice(0, 1), slice(None), slice(None)] if side == 0 else \
                 [slice(nx - 1, nx), slice(None), slice(None)]
        elif axis == 1:
            sl = [slice(None), slice(0, 1), slice(None)] if side == 0 else \
                 [slice(None), slice(ny - 1, ny), slice(None)]
        else:
            sl = [slice(None), slice(None), slice(0, 1)] if side == 0 else \
                 [slice(None), slice(None), slice(nz - 1, nz)]
        cells = idx[tuple(sl)].ravel()
        # Полуячейка: расстояние от центра ячейки до грани
        half = 0.5 * (grid.dx if axis == 0 else grid.dy if axis == 1 else grid.dz)
        area = (grid.dy * grid.dz) if axis == 0 else \
               (grid.dx * grid.dz) if axis == 1 else (grid.dx * grid.dy)
        k_face = k.ravel()[cells]
        t_ext = k_face * area / half
        np.add.at(diag, cells, t_ext)
        np.add.at(rhs, cells, t_ext * h_d)

    if not has_dirichlet:
        # Чистый Нейман: система вырождена по напору (произвольная константа)
        # и несогласована, если суммарный расход в блок не равен нулю.
        # Избыток/недостаток раствора равномерно распределяется по внешним
        # граням блока (модель «утечка через границы»), что делает задачу
        # согласованной и физически осмысленной.
        net = float(np.sum(rhs))
        meta["neumann_leakage"] = net
        if abs(net) > 1e-9:
            # вес ячейки по числу и площади выходящих на границу граней
            weight = np.zeros(n)
            for face, (axis, side) in FACE_AXES.items():
                if axis == 0:
                    sl = (0, slice(None), slice(None)) if side == 0 else (-1, slice(None), slice(None))
                elif axis == 1:
                    sl = (slice(None), 0, slice(None)) if side == 0 else (slice(None), -1, slice(None))
                else:
                    sl = (slice(None), slice(None), 0) if side == 0 else (slice(None), slice(None), -1)
                cells = idx[sl].ravel()
                np.add.at(weight, cells, 1.0)
            tot_w = weight.sum()
            if tot_w > 0:
                # уменьшаем источник в приграничных ячейках на величину избытка
                rhs = rhs - net * (weight / tot_w)

    rows.append(np.arange(n)); cols.append(np.arange(n)); data.append(diag)
    A = sp.coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n, n),
    ).tocsr()

    meta["has_dirichlet"] = has_dirichlet
    meta["diag"] = diag
    if not has_dirichlet:
        # Устранение вырожденности по константе: система A·x = b имеет
        # бесконечное семейство решений (x + const). Вводим условие
        # sum(x) = 0 через множитель Лагранжа:
        #     [ A   u ] [x]   [b]
        #     [ uᵀ  0 ] [λ] = [0],   u = 1/√n
        # Решение этой расширенной системы совпадает с решением
        # A·x = b при нулевом среднем (сдвиг напора не влияет на потоки).
        alpha = max(abs(diag).max(), 1.0) * 0.1
        u = sp.csr_matrix(np.ones((n, 1)) / np.sqrt(n))
        meta["A_ext"] = sp.bmat([[A, u], [u.T, None]], format="csc")
        meta["alpha"] = alpha
    return A, rhs, meta


def _face_perm(grid: GridSpec, k: np.ndarray, axis: int) -> np.ndarray:
    """Проницаемость на грани между соседними ячейками (гармоническое среднее)."""
    if axis == 0:
        return _harmonic(k[:-1, :, :], k[1:, :, :])
    if axis == 1:
        return _harmonic(k[:, :-1, :], k[:, 1:, :])
    return _harmonic(k[:, :, :-1], k[:, :, 1:])


def compute_fluxes(
    grid: GridSpec,
    k: np.ndarray,
    head: np.ndarray,
    kv_ratio: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Объёмные расходы Дарси (м2/сут) в центрах ячеек и объёмные расходы
    через грани (м3/сут).

    Расход через грань — величина, необходимая для адвективного переноса
    растворённых веществ: она должна вычисляться именно по грани, а не
    по осреднённому значению в ячейках.
    """
    nx, ny, nz = grid.shape
    qx = np.zeros(grid.shape)
    qy = np.zeros(grid.shape)
    qz = np.zeros(grid.shape)
    fx = np.zeros((nx - 1, ny, nz))
    fy = np.zeros((nx, ny - 1, nz))
    fz = np.zeros((nx, ny, nz - 1))

    for axis, (q, face, area, dist) in enumerate((
        (qx, fx, grid.dy * grid.dz, grid.dx),
        (qy, fy, grid.dx * grid.dz, grid.dy),
        (qz, fz, grid.dx * grid.dy, grid.dz),
    )):
        kf = _face_perm(grid, k, axis)
        if axis == 2:
            kf = kf * kv_ratio
        h = head
        if axis == 0:
            dh = (h[1:, :, :] - h[:-1, :, :]) / grid.dx
            fl = -kf * dh                       # расход Дарси, м/сут
            face[:, :, :] = fl * area           # объёмный расход через грань, м3/сут
            q[:-1, :, :] += 0.5 * fl
            q[1:, :, :] += 0.5 * fl
        elif axis == 1:
            dh = (h[:, 1:, :] - h[:, :-1, :]) / grid.dy
            fl = -kf * dh
            face[:, :, :] = fl * area
            q[:, :-1, :] += 0.5 * fl
            q[:, 1:, :] += 0.5 * fl
        else:
            dh = (h[:, :, 1:] - h[:, :, :-1]) / grid.dz
            fl = -kf * dh
            face[:, :, :] = fl * area
            q[:, :, :-1] += 0.5 * fl
            q[:, :, 1:] += 0.5 * fl
    return qx, qy, qz, fx, fy, fz


def solve(
    grid: GridSpec,
    k: np.ndarray,
    sources: np.ndarray,
    porosity: np.ndarray | None = None,
    bc_faces: dict[str, BoundaryCondition] | None = None,
    bc_heads: dict[str, float] | None = None,
    kv_ratio: float = 1.0,
) -> HydroResult:
    """Основной расчёт гидродинамики.

    sources — м3/сут в каждой ячейке (закачка +, отбор −).
    """
    k = np.asarray(k, dtype=float)
    # Защита от нулевых/отрицательных значений проницаемости
    k = np.where(np.isfinite(k) & (k > 0), k, 1e-6)

    A, rhs, meta = build_system(grid, k, sources, bc_faces, bc_heads, kv_ratio)
    try:
        if "A_ext" in meta:
            # Чистый Нейман: расширенная система с условием sum(x) = 0
            sol = spla.spsolve(meta["A_ext"], np.concatenate([rhs, [0.0]]))
            head = np.asarray(sol[:grid.n], dtype=float)
            solver = "spsolve+lagrange"
        else:
            head = spla.spsolve(A.tocsc(), rhs)
            solver = "spsolve"
        if not np.all(np.isfinite(head)):
            raise RuntimeError("non-finite")
    except Exception:
        head = spla.lsqr(A, rhs, atol=1e-10, btol=1e-10, iter_lim=5000)[0]
        solver = "lsqr"

    head = np.asarray(head, dtype=float).reshape(grid.shape)
    if not meta.get("has_dirichlet"):
        # Фиксация отсчёта напора: среднее равно нулю
        head = head - float(np.mean(head))
    qx, qy, qz, fx, fy, fz = compute_fluxes(grid, k, head, kv_ratio)

    phi = np.ones(grid.shape) if porosity is None else np.where(
        np.asarray(porosity) > 1e-6, np.asarray(porosity), 1e-6)
    vx, vy, vz = qx / phi, qy / phi, qz / phi
    speed = np.sqrt(vx**2 + vy**2 + vz**2)

    total_in = float(np.sum(sources[sources > 0]))
    total_out = float(-np.sum(sources[sources < 0]))

    # Баланс через граничные грани: расход, пересекающий внешние грани блока.
    # Положительный приток означает поступление воды в блок.
    b_in = 0.0
    b_out = 0.0
    for face, (axis, side) in FACE_AXES.items():
        cond_bc = (bc_faces or {}).get(face, BoundaryCondition.NEUMANN)
        if cond_bc != BoundaryCondition.DIRICHLET:
            continue
        h_d = float((bc_heads or {}).get(face, 0.0))
        if axis == 0:
            sl = 0 if side == 0 else -1
            area = grid.dy * grid.dz
            k_face = k[sl, :, :]
            h_face = head[sl, :, :]
            half = 0.5 * grid.dx
        elif axis == 1:
            sl = 0 if side == 0 else -1
            area = grid.dx * grid.dz
            k_face = k[:, sl, :]
            h_face = head[:, sl, :]
            half = 0.5 * grid.dy
        else:
            sl = 0 if side == 0 else -1
            area = grid.dx * grid.dy
            k_face = k[:, :, sl]
            h_face = head[:, :, sl]
            half = 0.5 * grid.dz
        k_eff = k_face * (kv_ratio if axis == 2 else 1.0)
        # приток извне внутрь блока (плюс — вода входит в блок)
        q_in = k_eff * (h_d - h_face) / half * area
        b_in += float(np.sum(q_in[q_in > 0]))
        b_out += float(-np.sum(q_in[q_in < 0]))

    residual = total_in - total_out - (b_in - b_out)

    return HydroResult(
        head=head, vx=vx, vy=vy, vz=vz, speed=speed,
        darcy_x=qx, darcy_y=qy, darcy_z=qz,
        face_x=fx, face_y=fy, face_z=fz,
        residual=residual,
        total_in=total_in, total_out=total_out,
        boundary_in=b_in, boundary_out=b_out,
        solver=solver,
    )


def build_sources(
    grid: GridSpec,
    wells,
    injection_factor: float = 1.0,
    production_factor: float = 1.0,
) -> tuple[np.ndarray, dict[str, list[tuple[int, int, int]]]]:
    """Распределение дебитов скважин по ячейкам фильтра.

    Возвращает источники (м3/сут) и карту «скважина → ячейки фильтра».
    """
    q = np.zeros(grid.shape)
    cells_map: dict[str, list[tuple[int, int, int]]] = {}
    for w in wells:
        cells = w.screen_cells(grid)
        cells_map[w.id] = cells
        if not cells or abs(w.rate) < 1e-12:
            continue
        from .models import WellType
        if w.well_type == WellType.INJECTION:
            total = w.rate * injection_factor
        elif w.well_type == WellType.PRODUCTION:
            total = -abs(w.rate) * production_factor
        else:
            continue
        for (i, j, k) in cells:
            q[i, j, k] += total / len(cells)
    return q, cells_map
