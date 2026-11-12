"""Базовые структуры данных проекта гидродинамического моделирования ПСВ урана.

Единицы измерения (SI-подобные, принятые в отрасли):
    координаты      — м (локальная система координат блока)
    глубина         — м (ось z направлена вниз, положительна вниз)
    проницаемость   — м/сут (коэффициент фильтрации k)
    пористость      — д.е.
    дебит           — м3/сут
    напор h         — м
    скорость        — м/сут
    кислотность     — г/л H2SO4
    уран в растворе — мг/л
    уран в руде     — кг/т (аналог), запас — кг
"""

from __future__ import annotations

import dataclasses
import datetime as _dt
import enum
from dataclasses import dataclass, field
from typing import Any

import numpy as np


class WellType(str, enum.Enum):
    INJECTION = "injection"   # закачка выщелачивающего раствора (ВР)
    PRODUCTION = "production" # откачка продуктивного раствора (ПР)
    OBSERVATION = "observation"


class BoundaryCondition(str, enum.Enum):
    NEUMANN = "neumann"  # непроницаемость
    DIRICHLET = "dirichlet"  # заданный напор на границе блока


@dataclass
class GridSpec:
    """Описание структурированной расчётной сетки (равномерной в каждом направлении)."""

    x0: float = 0.0
    y0: float = 0.0
    z0: float = 0.0
    nx: int = 40
    ny: int = 32
    nz: int = 12
    dx: float = 25.0
    dy: float = 25.0
    dz: float = 5.0

    @property
    def shape(self) -> tuple[int, int, int]:
        return (self.nx, self.ny, self.nz)

    @property
    def n_cells(self) -> int:
        return self.nx * self.ny * self.nz

    @property
    def cell_volume(self) -> float:
        return self.dx * self.dy * self.dz

    def centers(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Координаты центров ячеек: сетки (nx,ny,nz)."""
        x = self.x0 + (np.arange(self.nx) + 0.5) * self.dx
        y = self.y0 + (np.arange(self.ny) + 0.5) * self.dy
        z = self.z0 + (np.arange(self.nz) + 0.5) * self.dz
        return (
            np.broadcast_to(x[:, None, None], self.shape).copy(),
            np.broadcast_to(y[None, :, None], self.shape).copy(),
            np.broadcast_to(z[None, None, :], self.shape).copy(),
        )

    def bounds(self) -> dict[str, float]:
        return {
            "x_min": self.x0,
            "x_max": self.x0 + self.nx * self.dx,
            "y_min": self.y0,
            "y_max": self.y0 + self.ny * self.dy,
            "z_min": self.z0,
            "z_max": self.z0 + self.nz * self.dz,
        }


@dataclass
class ScreenInterval:
    """Интервал фильтрации скважины (глубины в м, z направлена вниз)."""

    top: float
    bottom: float

    def length(self) -> float:
        return max(0.0, self.bottom - self.top)


@dataclass
class TrajectoryPoint:
    x: float
    y: float
    z: float


@dataclass
class Well:
    """Скважина: траектория, интервал фильтра, тип и режим работы."""

    id: str
    name: str
    well_type: WellType = WellType.INJECTION
    trajectory: list[TrajectoryPoint] = field(default_factory=list)
    screen: ScreenInterval | None = None
    # Многоуровневые системы отработки: несколько интервалов фильтрации
    screens: list[ScreenInterval] = field(default_factory=list)
    # Дебит, м3/сут. Положительный — закачка (для закачивающих), отрицательный — отбор.
    rate: float = 0.0
    # Расход раствора, г/л H2SO4 на входе (для закачивающих скважин)
    acid_in: float = 20.0
    # Эксплуатация: интервал дат работы и изменение дебита
    schedule: list[tuple[str, float]] = field(default_factory=list)
    notes: str = ""

    def head_xy(self) -> tuple[float, float]:
        if self.trajectory:
            return self.trajectory[0].x, self.trajectory[0].y
        return 0.0, 0.0

    def screen_cells(self, grid: GridSpec) -> list[tuple[int, int, int]]:
        """Индексы ячеек, пересекающих интервалы фильтра (или ствол, если фильтр не задан)."""
        idx: list[tuple[int, int, int]] = []
        xc, yc, zc = grid.centers()
        if not self.trajectory:
            return idx
        intervals = list(self.screens) if self.screens else ([self.screen] if self.screen else [])
        if not intervals:
            return idx
        # Ствол аппроксимируется вертикальным отрезком в плане на средней отметке фильтра
        x_head, y_head = self.head_xy()
        sel_x = np.abs(xc - x_head) <= 0.5 * grid.dx + 1e-9
        sel_y = np.abs(yc - y_head) <= 0.5 * grid.dy + 1e-9
        for iv in intervals:
            # все ячейки, попадающие в интервал фильтрации
            sel_z = (zc >= iv.top - 1e-9) & (zc <= iv.bottom + 1e-9)
            mask = sel_z & sel_x & sel_y
            for i, j, k in np.argwhere(mask):
                idx.append((int(i), int(j), int(k)))
        return idx


@dataclass
class RateStage:
    """Режим (этап) работы: интервал дат и заданная нагрузка на скважины."""

    name: str
    start: str
    end: str
    # Множители дебитов относительно базового расхода скважин
    injection_factor: float = 1.0
    production_factor: float = 1.0

    def days(self) -> float:
        d0 = _dt.date.fromisoformat(self.start)
        d1 = _dt.date.fromisoformat(self.end)
        return max(1.0, (d1 - d0).days)


@dataclass
class ChemicalParams:
    """Параметры химической кинетики (массопереноса)."""

    acid_concentration: float = 20.0       # г/л H2SO4 (реагент по умолчанию)
    use_acid_data: bool = True            # учитывать загруженные данные по кислотности
    fe3_conc: float = 0.0                 # г/л Fe(III)
    feoh3_conc: float = 0.35              # масс. % Fe(OH)3
    caco3_conc: float = 0.0               # % CaCO3
    integration_days: float = 365.0       # время интегрирования, сут
    k1: float = 0.0003                    # скорость кислотного выщелачивания UO2, 1/сут
    k2: float = 0.00005                   # скорость восстановления Fe(III), 1/сут
    k3: float = 0.00002                   # осаждение урана (потери), 1/сут
    k4: float = 0.00005                   # сорбция урана на глине, 1/сут
    beta: float = 0.0                     # коэффициент замедления из-за кольматации
    flowline_k1: float | None = None
    flowline_k2: float | None = None
    flowline_k3: float | None = None
    flowline_k4: float | None = None
    flowline_beta: float | None = None
    flowline_acid: float = 20.0
    # Температурная поправка (Аррениус) для скоростей реакций
    temperature: float = 20.0
    activation_energy: float = 20.0

    def k1_fl(self) -> float:
        return self.flowline_k1 if self.flowline_k1 is not None else self.k1

    def k2_fl(self) -> float:
        return self.flowline_k2 if self.flowline_k2 is not None else self.k2

    def k3_fl(self) -> float:
        return self.flowline_k3 if self.flowline_k3 is not None else self.k3

    def k4_fl(self) -> float:
        return self.flowline_k4 if self.flowline_k4 is not None else self.k4


@dataclass
class PropertyModel:
    """Свойства пласта, заданные по скважинам и интерполируемые на сетку."""

    permeability: list[tuple[float, float, float, float]] = field(default_factory=list)
    # (x, y, z, k) — коэффициент фильтрации, м/сут
    porosity: list[tuple[float, float, float, float]] = field(default_factory=list)
    uranium: list[tuple[float, float, float, float]] = field(default_factory=list)
    # (x, y, z, содержание урана, кг/т)
    clay: list[tuple[float, float, float, float]] = field(default_factory=list)
    # (x, y, z, доля глинистой компоненты, 0..1)
    lithology: list[tuple[float, float, float, str]] = field(default_factory=list)
    # (x, y, z, код литологии)

    def all_points(self) -> list[tuple[float, float, float]]:
        pts = [p[:3] for p in self.permeability]
        pts += [p[:3] for p in self.uranium]
        return pts


@dataclass
class Project:
    """Проект моделирования: сетка, скважины, свойства, результаты."""

    name: str = "Новый проект"
    created: str = field(default_factory=lambda: _dt.datetime.now().isoformat(timespec="seconds"))
    grid: GridSpec = field(default_factory=GridSpec)
    wells: list[Well] = field(default_factory=list)
    props: PropertyModel = field(default_factory=PropertyModel)
    chem: ChemicalParams = field(default_factory=ChemicalParams)
    boundary: BoundaryCondition = BoundaryCondition.NEUMANN
    boundary_head: float = 0.0
    # Условия по каждой грани блока: "neumann" | "dirichlet"
    bc_faces: dict[str, str] = field(default_factory=dict)
    bc_heads: dict[str, float] = field(default_factory=dict)
    stages: list[RateStage] = field(default_factory=list)
    # Расчётные поля (numpy-массивы), заполняются после расчёта
    results: dict[str, Any] = field(default_factory=dict)
    # Расчётные поля в виде numpy-массивов (не сериализуются в JSON)
    arrays: dict[str, Any] = field(default_factory=dict)

    def well_by_id(self, wid: str) -> Well | None:
        for w in self.wells:
            if w.id == wid:
                return w
        return None

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "created": self.created,
            "grid": dataclasses.asdict(self.grid),
            "boundary": self.boundary.value,
            "boundary_head": self.boundary_head,
            "bc_faces": self.bc_faces,
            "bc_heads": self.bc_heads,
            "chem": dataclasses.asdict(self.chem),
            "stages": [dataclasses.asdict(s) for s in self.stages],
            "wells": [
                {
                    "id": w.id,
                    "name": w.name,
                    "well_type": w.well_type.value,
                    "trajectory": [dataclasses.asdict(p) for p in w.trajectory],
                    "screen": dataclasses.asdict(w.screen) if w.screen else None,
                    "screens": [dataclasses.asdict(s) for s in w.screens],
                    "rate": w.rate,
                    "acid_in": w.acid_in,
                    "schedule": w.schedule,
                    "notes": w.notes,
                }
                for w in self.wells
            ],
            "props": {
                "permeability": self.props.permeability,
                "porosity": self.props.porosity,
                "uranium": self.props.uranium,
                "clay": self.props.clay,
                "lithology": self.props.lithology,
            },
            "results": {
                k: (v.tolist() if isinstance(v, np.ndarray) else v)
                for k, v in self.results.items()
                if not k.startswith("_")
            },
        }
