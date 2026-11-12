"""Генерация демонстрационного проекта.

Синтетический, но физически осмысленный блок, повторяющий типовое
размещение скважин при подземном скважинном выщелачивании: чередующиеся
ряды закачных и откачных скважин, две отработки (верхняя и нижняя),
неоднородный коэффициент фильтрации, содержания урана и литология.

Все данные синтетические и предназначены исключительно для проверки
работоспособности программы.
"""

from __future__ import annotations

import numpy as np

from .models import (
    BoundaryCondition,
    ChemicalParams,
    GridSpec,
    ScreenInterval,
    TrajectoryPoint,
    Well,
    WellType,
)

BLOCK_X = 800.0
BLOCK_Y = 600.0
BLOCK_Z = 60.0
CELL = 25.0


def build_demo_project(project) -> None:
    rng = np.random.default_rng(20240115)

    project.grid = GridSpec(x0=0, y0=0, z0=0, nx=32, ny=24, nz=6,
                            dx=CELL, dy=CELL, dz=10.0)
    project.boundary = BoundaryCondition.NEUMANN
    project.bc_faces = {f: "neumann" for f in
                        ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")}
    project.bc_heads = {}
    project.chem = ChemicalParams(
        acid_concentration=20.0, use_acid_data=True,
        fe3_conc=2.0, feoh3_conc=0.35, caco3_conc=0.0,
        integration_days=365.0,
        k1=0.0003, k2=0.00005, k3=0.00002, k4=0.00005,
        beta=0.6, temperature=22.0, activation_energy=20.0,
        flowline_k1=0.0003, flowline_k2=0.00005, flowline_k3=0.00002,
        flowline_k4=0.00005, flowline_beta=0.6, flowline_acid=20.0,
    )

    # ---- сеть скважин: ряды закачки/откачки ----
    wells: list[Well] = []
    n_rows = 6          # чередование рядов закачки/откачки
    per_row = 7
    x_start, x_step = 100.0, 100.0
    y_start, y_step = 100.0, 90.0

    idx = 0
    for row in range(n_rows):
        is_injection = (row % 2 == 0)
        for c in range(per_row):
            x = x_start + c * x_step
            y = y_start + row * y_step
            z_head = 2.0 + rng.uniform(-1.5, 1.5)
            idx += 1
            wid = f"W{idx:02d}"
            traj = [TrajectoryPoint(x=x, y=y, z=z_head)]
            # Нижняя отработка: 30-45 м; верхняя: 12-24 м
            screens = [ScreenInterval(top=30.0, bottom=45.0)]
            if row in (0, 5):
                screens = [ScreenInterval(top=12.0, bottom=24.0),
                           ScreenInterval(top=30.0, bottom=45.0)]
            rate = 620.0 if is_injection else -600.0
            rate *= float(rng.uniform(0.9, 1.1))
            w = Well(
                id=wid,
                name=f"{'З' if is_injection else 'О'}-{idx:02d}",
                well_type=WellType.INJECTION if is_injection else WellType.PRODUCTION,
                trajectory=traj,
                screens=screens,
                screen=screens[0],
                rate=rate,
                acid_in=20.0 if is_injection else 0.0,
                notes="Учебные данные",
            )
            wells.append(w)
    project.wells = wells

    # Баланс расходов: суммарный отбор приравнивается к суммарной закачке,
    # иначе при замкнутых границах задача гидродинамики несогласована.
    inj_total = sum(w.rate for w in wells if w.well_type == WellType.INJECTION)
    prod = [w for w in wells if w.well_type == WellType.PRODUCTION]
    if prod and inj_total > 0:
        for w in prod:
            w.rate = -inj_total / len(prod)

    # ---- свойства по скважинам ----
    k_all: list[tuple[float, float, float, float]] = []
    u_all: list[tuple[float, float, float, float]] = []
    p_all: list[tuple[float, float, float, float]] = []
    c_all: list[tuple[float, float, float, float]] = []
    l_all: list[tuple[float, float, float, str]] = []

    # Две «геологические» зоны с разной проницаемостью (тектонический разлом)
    def zone_k(x: float, y: float) -> float:
        base = 1.6
        if x > 480.0:
            base *= 0.45          # слабопроницаемый блок
        if 250.0 < y < 330.0:
            base *= 1.8           # русло, повышенная проницаемость
        return base

    for w in wells:
        x, y = w.head_xy()
        for s in w.screens:
            z_mid = 0.5 * (s.top + s.bottom)
            k_val = zone_k(x, y) * float(rng.uniform(0.75, 1.3))
            k_all.append((x, y, z_mid, k_val))
            # Содержание урана: рудные «линзы»
            g = 0.55 * float(rng.uniform(0.4, 1.6))
            if 300.0 < x < 600.0 and 150.0 < y < 450.0:
                g *= 1.7
            u_all.append((x, y, z_mid, g))
            p_all.append((x, y, z_mid, float(rng.uniform(0.13, 0.26))))
            clay = 0.18 if x > 480.0 else 0.07
            c_all.append((x, y, z_mid, clay))
            lith = "глинистый алевролит" if x > 480.0 else "песчаник"
            l_all.append((x, y, z_mid, lith))

    # Дополнительные точки для более равномерного покрытия блока
    for _ in range(40):
        x = float(rng.uniform(30, BLOCK_X - 30))
        y = float(rng.uniform(30, BLOCK_Y - 30))
        z = float(rng.choice([18.0, 37.0]))
        k_all.append((x, y, z, zone_k(x, y) * float(rng.uniform(0.8, 1.25))))
        u_all.append((x, y, z, 0.5 * float(rng.uniform(0.3, 1.5))))
        p_all.append((x, y, z, float(rng.uniform(0.13, 0.25))))
        c_all.append((x, y, z, 0.15 if x > 480 else 0.07))
        l_all.append((x, y, z, "глинистый алевролит" if x > 480 else "песчаник"))

    project.props.permeability = k_all
    project.props.uranium = u_all
    project.props.porosity = p_all
    project.props.clay = c_all
    project.props.lithology = l_all
