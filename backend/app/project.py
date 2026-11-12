"""Оркестратор расчётного конвейера проекта.

Связывает модули: геология (интерполяция) → гидродинамика → линии тока →
химическая кинетика (массоперенос) → подсчёт запасов.
"""

from __future__ import annotations

import json
import os
import shutil
import uuid
import zipfile
from dataclasses import dataclass
from typing import Any

import numpy as np

from . import chemistry, flowlines, geo_ops, hydro, reserves
from .geostat import GeoSettings, Variogram, interpolate, nearest_lithology
from .models import (
    BoundaryCondition,
    ChemicalParams,
    GridSpec,
    Project,
    PropertyModel,
    TrajectoryPoint,
    Well,
    WellType,
)


@dataclass
class RunOptions:
    run_hydro: bool = True
    run_flowlines: bool = True
    run_chemistry: bool = True
    run_reserves: bool = True
    interpolation: str = "kriging"          # kriging | idw
    variogram_model: str = "spherical"
    variogram_sill: float = 1.0
    variogram_range: float = 150.0
    variogram_nugget: float = 0.0
    anisotropy: tuple[float, float, float] = (1.0, 1.0, 1.0)
    kv_ratio: float = 0.2                # анизотропия проницаемости по вертикали (k_z/k_xy)
    stagnation_threshold: float = 0.05       # м/сут (ниже — застойная зона)
    max_streamlines: int = 400
    contour_polygon: list[list[float]] | None = None
    acid_threshold: float = 5.0
    output_every: float = 30.0              # сут
    reactions: tuple[bool, bool, bool, bool] = (True, True, True, True)


def _geo_settings(opts: RunOptions) -> GeoSettings:
    return GeoSettings(
        method=opts.interpolation,
        variogram=Variogram(
            model=opts.variogram_model,
            sill=opts.variogram_sill,
            range_=opts.variogram_range,
            nugget=opts.variogram_nugget,
        ),
        anisotropy=opts.anisotropy,
    )


def build_property_models(project: Project, opts: RunOptions) -> dict[str, np.ndarray]:
    """Интерполяция скважинных свойств на расчётную сетку."""
    grid = project.grid
    xc, yc, zc = grid.centers()
    targets = np.stack([xc.ravel(), yc.ravel(), zc.ravel()], axis=1)
    settings = _geo_settings(opts)

    def interp(points, values, default, clip_to_data=True):
        if not points:
            return np.full(grid.n_cells, default)
        arr = np.array([p[3] for p in points], dtype=float)
        pts = np.array([p[:3] for p in points], dtype=float)
        res = interpolate(pts, arr, targets, settings)
        res = np.where(np.isfinite(res), res, default)
        if clip_to_data and arr.size:
            # Отсечение выбросов интерполяции диапазоном исходных данных:
            # физически проницаемость и содержание не могут выходить за
            # пределы, наблюдаемые в скважинах.
            lo = float(np.nanpercentile(arr, 1))
            hi = float(np.nanpercentile(arr, 99))
            if hi <= lo:
                lo, hi = float(np.min(arr)), float(np.max(arr))
            res = np.clip(res, lo, hi)
        return res.reshape(grid.shape)

    k = interp(project.props.permeability, None, 1.0)
    porosity = interp(project.props.porosity, None, 0.15)
    uranium = interp(project.props.uranium, None, 0.0)
    clay = interp(project.props.clay, None, 0.1)

    lith = np.full(grid.n_cells, "", dtype=object)
    if project.props.lithology:
        pts = np.array([p[:3] for p in project.props.lithology], dtype=float)
        labels = [p[3] for p in project.props.lithology]
        lith = nearest_lithology(pts, labels, targets)

    return {
        "k": k,
        "porosity": porosity,
        "uranium": uranium,
        "clay": clay,
        "lithology": lith.reshape(grid.shape),
    }


def _bc_maps(project: Project) -> tuple[dict[str, BoundaryCondition], dict[str, float]]:
    faces: dict[str, BoundaryCondition] = {}
    heads: dict[str, float] = {}
    for face in hydro.FACE_AXES:
        val = project.bc_faces.get(face)
        if val:
            faces[face] = BoundaryCondition(val)
        elif project.boundary == BoundaryCondition.DIRICHLET:
            faces[face] = BoundaryCondition.DIRICHLET
        if face in faces and faces[face] == BoundaryCondition.DIRICHLET:
            heads[face] = float(project.bc_heads.get(face, project.boundary_head))
    return faces, heads


def run_pipeline(project: Project, opts: RunOptions | None = None) -> dict[str, Any]:
    """Полный расчёт по проекту. Результаты сохраняются в project.results."""
    opts = opts or RunOptions()
    grid = project.grid
    out: dict[str, Any] = {}

    # 1. Геологические модели
    props = build_property_models(project, opts)
    out["property_stats"] = {
        "k_min": float(np.min(props["k"])), "k_max": float(np.max(props["k"])),
        "k_mean": float(np.mean(props["k"])),
        "poro_mean": float(np.mean(props["porosity"])),
        "u_max": float(np.max(props["uranium"])),
        "u_mean": float(np.mean(props["uranium"])),
        "points": len(project.props.permeability),
    }

    # 2. Гидродинамика
    if opts.run_hydro:
        sources, cells_map = hydro.build_sources(grid, project.wells)
        faces, heads = _bc_maps(project)
        hr = hydro.solve(grid, props["k"], sources, props["porosity"], faces, heads,
                         kv_ratio=opts.kv_ratio)
        out["hydro"] = {
            "residual": hr.residual,
            "total_in": hr.total_in,
            "total_out": hr.total_out,
            "boundary_in": hr.boundary_in,
            "boundary_out": hr.boundary_out,
            "head_min": float(np.min(hr.head)),
            "head_max": float(np.max(hr.head)),
            "speed_max": float(np.max(hr.speed)),
            "solver": hr.solver,
            "balance_ok": abs(hr.residual) < max(1e-6, 0.01 * max(hr.total_in, hr.total_out, 1.0)),
        }
        out["_hydro_result"] = hr
        out["_cells_map"] = cells_map

        # 3. Линии тока
        if opts.run_flowlines:
            seeds, labels = flowlines.seeds_from_wells(project.wells, grid)
            stops = flowlines.stop_cells_from_wells(project.wells, grid)
            if len(seeds) > opts.max_streamlines:
                sel = np.linspace(0, len(seeds) - 1, opts.max_streamlines).astype(int)
                seeds, labels = seeds[sel], [labels[i] for i in sel]
            sl = flowlines.trace(
                grid, hr.vx, hr.vy, hr.vz, seeds, labels, stops,
                stagnation_threshold=opts.stagnation_threshold,
                max_residence=max(365.0, project.chem.integration_days * 3),
            )
            out["streamlines"] = [s.to_dict() for s in sl]
            out["_streamlines"] = sl
            out["stagnant_fraction"] = float(np.mean(hr.speed < opts.stagnation_threshold))

    # 4. Массоперенос
    if opts.run_chemistry and "_hydro_result" in out:
        hr = out["_hydro_result"]
        times = list(np.arange(opts.output_every, project.chem.integration_days + 1e-9, opts.output_every))
        if not times or times[-1] < project.chem.integration_days:
            times.append(project.chem.integration_days)
        mt = chemistry.simulate_3d(
            grid, props["k"], props["porosity"], props["uranium"], props["clay"],
            hr.face_x, hr.face_y, hr.face_z,
            project.wells, project.chem, times, reactions=opts.reactions,
        )
        out["_mass_transfer"] = mt
        out["curves"] = {
            "time": mt.times,
            "pr_uranium": mt.pr_uranium,
            "vr_acid": mt.vr_acid,
            "cumulative_u": mt.cumulative_u,
            "cumulative_acid": mt.cumulative_acid,
            "prod_volume": mt.prod_volume,
            "inj_volume": mt.inj_volume,
        }

        # 5. Запасы
        if opts.run_reserves:
            contour = None
            if opts.contour_polygon:
                contour = reserves.contour_from_polygon(grid, opts.contour_polygon)
            rr = reserves.calculate(
                grid, mt.residual[-1], mt.acid[-1],
                acid_threshold=opts.acid_threshold,
                contour_mask=contour,
                porosity=props["porosity"],
            )
            out["reserves"] = {
                "reserve_kg": rr.reserve_kg, "reserve_t": rr.reserve_t,
                "acidified_volume": rr.acidified_volume, "acidified_cells": rr.acidified_cells,
                "area_m2": rr.area_m2, "area_ha": rr.area_ha,
                "avg_grade": rr.avg_grade, "avg_grade_pct": rr.avg_grade_pct,
                "avg_thickness": rr.avg_thickness,
                "weighted_metre_percent": rr.weighted_metre_percent,
                "min_grade": rr.min_grade, "max_grade": rr.max_grade,
                "total_ore_t": rr.total_ore_t,
                "contour_area_m2": rr.contour_area_m2,
                "contour_volume": rr.contour_volume,
                "by_layer": rr.by_layer,
            }

    # Очищаем несериализуемые объекты
    clean = {k: v for k, v in out.items() if not k.startswith("_")}
    project.results = dict(clean)          # копия: служебные ключи не попадут в ответ API
    project.arrays = {
        "k": props["k"], "porosity": props["porosity"],
        "uranium": props["uranium"], "clay": props["clay"],
    }
    if "_hydro_result" in out:
        hr = out["_hydro_result"]
        project.arrays.update({
            "head": hr.head, "speed": hr.speed,
            "vx": hr.vx, "vy": hr.vy, "vz": hr.vz,
        })
    if "_mass_transfer" in out:
        mt = out["_mass_transfer"]
        project.arrays.update({
            "acid_final": mt.acid[-1], "uranium_final": mt.uranium[-1],
            "residual_final": mt.residual[-1],
        })
        project.results["_times"] = mt.times
        project.results["_acid_history"] = mt.acid
        project.results["_uranium_history"] = mt.uranium
    if "_streamlines" in out:
        project.results["_streamlines"] = out["_streamlines"]
    return clean


# ---------------------------------------------------------------- хранение

STORE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "data", "projects")


class ProjectStore:
    """Хранилище проектов: в памяти + сохранение на диск."""

    def __init__(self, base_dir: str | None = None):
        self.base_dir = base_dir or STORE_DIR
        os.makedirs(self.base_dir, exist_ok=True)
        self._projects: dict[str, Project] = {}
        self._options: dict[str, RunOptions] = {}

    def create(self, name: str = "Новый проект") -> Project:
        pid = uuid.uuid4().hex[:12]
        p = Project(name=name)
        self._projects[pid] = p
        self._options[pid] = RunOptions()
        return p

    def get(self, pid: str) -> Project | None:
        return self._projects.get(pid)

    def list(self) -> list[dict]:
        out = []
        for pid, p in self._projects.items():
            out.append({
                "id": pid, "name": p.name, "created": p.created,
                "wells": len(p.wells),
                "grid": p.grid.shape,
            })
        return out

    def delete(self, pid: str) -> bool:
        self._projects.pop(pid, None)
        self._options.pop(pid, None)
        path = os.path.join(self.base_dir, f"{pid}.s3dp")
        if os.path.exists(path):
            os.remove(path)
        return True

    def options(self, pid: str) -> RunOptions:
        return self._options.setdefault(pid, RunOptions())

    def save(self, pid: str) -> str:
        p = self._projects[pid]
        path = os.path.join(self.base_dir, f"{pid}.s3dp")
        arrays: dict[str, np.ndarray] = {}
        for k, v in p.arrays.items():
            arrays[k] = np.asarray(v)
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("project.json", json.dumps(p.to_dict(), ensure_ascii=False, indent=1))
            if arrays:
                import io
                buf = io.BytesIO()
                np.savez_compressed(buf, **arrays)
                z.writestr("arrays.npz", buf.getvalue())
        return path

    def load(self, path: str) -> str:
        with zipfile.ZipFile(path, "r") as z:
            data = json.loads(z.read("project.json").decode("utf-8"))
            p = Project(name=data["name"], created=data.get("created"))
            g = data["grid"]
            p.grid = GridSpec(**g)
            p.boundary = BoundaryCondition(data.get("boundary", "neumann"))
            p.boundary_head = float(data.get("boundary_head", 0.0))
            p.bc_faces = data.get("bc_faces", {}) or {}
            p.bc_heads = data.get("bc_heads", {}) or {}
            p.chem = ChemicalParams(**data.get("chem", {}))
            pr = data.get("props", {}) or {}
            p.props = PropertyModel(
                permeability=[tuple(r) for r in pr.get("permeability", [])],
                porosity=[tuple(r) for r in pr.get("porosity", [])],
                uranium=[tuple(r) for r in pr.get("uranium", [])],
                clay=[tuple(r) for r in pr.get("clay", [])],
                lithology=[tuple(r) for r in pr.get("lithology", [])],
            )
            for w in data.get("wells", []):
                traj = [TrajectoryPoint(**t) for t in w.get("trajectory", [])]
                screen = None
                if w.get("screen"):
                    screen = geo_ops.make_screen(w["screen"])
                screens = [geo_ops.make_screen(s) for s in w.get("screens", []) or []]
                p.wells.append(Well(
                    id=w["id"], name=w["name"],
                    well_type=WellType(w.get("well_type", "injection")),
                    trajectory=traj, screen=screen, screens=screens,
                    rate=float(w.get("rate", 0.0)), acid_in=float(w.get("acid_in", 20.0)),
                    schedule=[tuple(s) for s in w.get("schedule", [])],
                    notes=w.get("notes", ""),
                ))
            if "arrays.npz" in z.namelist():
                import io
                buf = io.BytesIO(z.read("arrays.npz"))
                with np.load(buf) as npz:
                    p.arrays = {k: npz[k] for k in npz.files}
        pid = uuid.uuid4().hex[:12]
        self._projects[pid] = p
        self._options[pid] = RunOptions()
        return pid
