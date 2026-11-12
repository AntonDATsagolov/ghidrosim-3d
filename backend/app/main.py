"""REST API и веб-сервер программы «ГидроСимулятор 3D»."""

from __future__ import annotations

import io
import json
import os
import time
import uuid

import numpy as np
from fastapi import Body, FastAPI, File, HTTPException, Query, UploadFile
from fastapi.responses import FileResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from . import importers, reserves
from .models import (
    BoundaryCondition,
    ChemicalParams,
    GridSpec,
    Project,
    PropertyModel,
    ScreenInterval,
    TrajectoryPoint,
    Well,
    WellType,
)
from .project import ProjectStore, RunOptions, run_pipeline

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FRONTEND_DIR = os.path.join(BASE_DIR, "frontend")

app = FastAPI(title="ГидроСимулятор 3D", version="1.0.0")
store = ProjectStore()


# ------------------------------------------------------------------ схемы

class GridIn(BaseModel):
    x0: float = 0.0
    y0: float = 0.0
    z0: float = 0.0
    nx: int = Field(40, ge=2, le=400)
    ny: int = Field(32, ge=2, le=400)
    nz: int = Field(12, ge=1, le=200)
    dx: float = Field(25.0, gt=0)
    dy: float = Field(25.0, gt=0)
    dz: float = Field(5.0, gt=0)


class WellIn(BaseModel):
    id: str | None = None
    name: str
    well_type: str = "injection"
    x: float
    y: float
    z: float = 0.0
    screen_top: float | None = None
    screen_bottom: float | None = None
    rate: float = 0.0
    acid_in: float = 20.0
    trajectory: list[dict] | None = None
    notes: str = ""


class PropsIn(BaseModel):
    permeability: list[list[float]] = Field(default_factory=list)
    porosity: list[list[float]] = Field(default_factory=list)
    uranium: list[list[float]] = Field(default_factory=list)
    clay: list[list[float]] = Field(default_factory=list)
    lithology: list[list] = Field(default_factory=list)


class BCIn(BaseModel):
    boundary: str = "neumann"
    boundary_head: float = 0.0
    bc_faces: dict[str, str] = Field(default_factory=dict)
    bc_heads: dict[str, float] = Field(default_factory=dict)


class RunIn(BaseModel):
    run_hydro: bool = True
    run_flowlines: bool = True
    run_chemistry: bool = True
    run_reserves: bool = True
    interpolation: str = "kriging"
    variogram_model: str = "spherical"
    variogram_sill: float = 1.0
    variogram_range: float = 150.0
    variogram_nugget: float = 0.0
    anisotropy: list[float] = Field(default_factory=lambda: [1.0, 1.0, 1.0])
    kv_ratio: float = 0.2
    stagnation_threshold: float = 0.05
    max_streamlines: int = 400
    contour_polygon: list[list[float]] | None = None
    acid_threshold: float = 5.0
    output_every: float = 30.0
    reactions: list[bool] = Field(default_factory=lambda: [True, True, True, True])


class ImportApplyIn(BaseModel):
    rows: list[dict]
    mapping: dict[str, str | None] = Field(default_factory=dict)
    as_properties: bool = True


# ------------------------------------------------------------------ утилиты

def _well_from_in(w: WellIn) -> Well:
    wid = w.id or uuid.uuid4().hex[:8]
    traj = [TrajectoryPoint(**t) for t in (w.trajectory or [])]
    if not traj:
        traj = [TrajectoryPoint(x=w.x, y=w.y, z=w.z)]
    screen = None
    if w.screen_top is not None and w.screen_bottom is not None and w.screen_bottom > w.screen_top:
        screen = ScreenInterval(top=w.screen_top, bottom=w.screen_bottom)
    return Well(
        id=wid, name=w.name,
        well_type=WellType(w.well_type) if w.well_type in {t.value for t in WellType} else WellType.INJECTION,
        trajectory=traj, screen=screen, rate=w.rate, acid_in=w.acid_in, notes=w.notes,
    )


def _project_summary(p: Project) -> dict:
    return {
        "id": None,
        "name": p.name,
        "created": p.created,
        "grid": {
            **p.grid.bounds(),
            "x0": p.grid.x0, "y0": p.grid.y0, "z0": p.grid.z0,
            "nx": p.grid.nx, "ny": p.grid.ny, "nz": p.grid.nz,
            "dx": p.grid.dx, "dy": p.grid.dy, "dz": p.grid.dz,
            "shape": p.grid.shape,
        },
        "wells": [
            {
                "id": w.id, "name": w.name, "well_type": w.well_type.value,
                "x": w.head_xy()[0], "y": w.head_xy()[1],
                "z": w.trajectory[0].z if w.trajectory else 0.0,
                "screen_top": w.screen.top if w.screen else None,
                "screen_bottom": w.screen.bottom if w.screen else None,
                "screens": [{"top": s.top, "bottom": s.bottom} for s in w.screens],
                "rate": w.rate, "acid_in": w.acid_in, "notes": w.notes,
            }
            for w in p.wells
        ],
        "boundary": p.boundary.value,
        "bc_faces": p.bc_faces,
        "bc_heads": p.bc_heads,
        "chem": p.chem.__dict__,
        "property_counts": {
            "permeability": len(p.props.permeability),
            "uranium": len(p.props.uranium),
            "porosity": len(p.props.porosity),
            "clay": len(p.props.clay),
            "lithology": len(p.props.lithology),
        },
    }


def _need(pid: str) -> Project:
    p = store.get(pid)
    if p is None:
        raise HTTPException(404, "Проект не найден")
    return p


# ------------------------------------------------------------------ маршруты

@app.get("/api/health")
def health():
    return {"ok": True, "version": app.version}


@app.post("/api/projects")
def create_project(payload: dict = Body(default={})):
    p = store.create(payload.get("name", "Новый проект"))
    if payload.get("grid"):
        g = payload["grid"]
        p.grid = GridSpec(
            x0=g.get("x0", 0), y0=g.get("y0", 0), z0=g.get("z0", 0),
            nx=int(g.get("nx", 40)), ny=int(g.get("ny", 32)), nz=int(g.get("nz", 12)),
            dx=float(g.get("dx", 25)), dy=float(g.get("dy", 25)), dz=float(g.get("dz", 5)),
        )
    s = _project_summary(p)
    s["id"] = _pid_of(p)
    return s


def _pid_of(p: Project) -> str:
    for pid, pr in store._projects.items():
        if pr is p:
            return pid
    return ""


@app.get("/api/projects")
def list_projects():
    return {"projects": store.list()}


@app.get("/api/projects/{pid}")
def get_project(pid: str):
    p = _need(pid)
    s = _project_summary(p)
    s["id"] = pid
    if p.results:
        s["results_available"] = True
        s["hydro"] = p.results.get("hydro")
        s["reserves"] = p.results.get("reserves")
        s["curves"] = p.results.get("curves")
        s["property_stats"] = p.results.get("property_stats")
        s["stagnant_fraction"] = p.results.get("stagnant_fraction")
    return s


@app.delete("/api/projects/{pid}")
def delete_project(pid: str):
    _need(pid)
    store.delete(pid)
    return {"ok": True}


@app.put("/api/projects/{pid}/grid")
def set_grid(pid: str, g: GridIn):
    p = _need(pid)
    p.grid = GridSpec(**g.model_dump())
    return _project_summary(p)


@app.put("/api/projects/{pid}/wells")
def set_wells(pid: str, wells: list[WellIn]):
    p = _need(pid)
    p.wells = [_well_from_in(w) for w in wells]
    return _project_summary(p)


@app.put("/api/projects/{pid}/properties")
def set_properties(pid: str, props: PropsIn):
    p = _need(pid)
    p.props = PropertyModel(
        permeability=[tuple(r) for r in props.permeability],
        porosity=[tuple(r) for r in props.porosity],
        uranium=[tuple(r) for r in props.uranium],
        clay=[tuple(r) for r in props.clay],
        lithology=[tuple(r) for r in props.lithology],
    )
    return _project_summary(p)


@app.put("/api/projects/{pid}/chemistry")
def set_chemistry(pid: str, chem: ChemicalParams):
    p = _need(pid)
    p.chem = chem
    return {"ok": True}


@app.put("/api/projects/{pid}/bc")
def set_bc(pid: str, bc: BCIn):
    p = _need(pid)
    p.boundary = BoundaryCondition(bc.boundary)
    p.boundary_head = bc.boundary_head
    p.bc_faces = bc.bc_faces
    p.bc_heads = bc.bc_heads
    return {"ok": True}


@app.post("/api/projects/{pid}/import")
async def import_file(pid: str, file: UploadFile = File(...)):
    _need(pid)
    data = await file.read()
    df = importers.read_table(data, file.filename or "data.csv")
    mapping = importers.auto_map_columns(df)
    rows = importers.rows_to_wells(df, mapping)
    return {
        "filename": file.filename,
        "columns": list(df.columns),
        "mapping": mapping,
        "rows": rows,
        "count": len(rows),
        "preview": df.head(10).astype(object).where(pd_notna(df), None).to_dict(orient="records")
        if False else json.loads(df.head(10).to_json(orient="records", force_ascii=False)),
    }


def pd_notna(df):
    import pandas as pd
    return pd.notna(df)


@app.post("/api/projects/{pid}/import/apply")
def import_apply(pid: str, payload: ImportApplyIn):
    p = _need(pid)
    mapping = payload.mapping or {}
    # Новые скважины
    for r in payload.rows:
        wi = WellIn(
            name=str(r.get("name") or "скважина"),
            well_type=r.get("well_type") or "injection",
            x=float(r.get("x") or 0.0), y=float(r.get("y") or 0.0),
            z=float(r.get("z") or 0.0),
            screen_top=r.get("screen_top"), screen_bottom=r.get("screen_bottom"),
            rate=float(r.get("rate") or 0.0),
            acid_in=float(r.get("acid_in") or 20.0),
        )
        p.wells.append(_well_from_in(wi))
    # Свойства пласта
    if payload.as_properties:
        for r in payload.rows:
            x, y, z = r.get("x"), r.get("y"), r.get("z")
            if x is None or y is None:
                continue
            if r.get("k") is not None:
                p.props.permeability.append((float(x), float(y), float(z), float(r["k"])))
            if r.get("uranium") is not None:
                p.props.uranium.append((float(x), float(y), float(z), float(r["uranium"])))
            if r.get("porosity") is not None:
                p.props.porosity.append((float(x), float(y), float(z), float(r["porosity"])))
            if r.get("clay") is not None:
                p.props.clay.append((float(x), float(y), float(z), float(r["clay"])))
            if r.get("lithology"):
                p.props.lithology.append((float(x), float(y), float(z), str(r["lithology"])))
    return _project_summary(p)


@app.post("/api/projects/{pid}/run")
def run(pid: str, opts: RunIn):
    p = _need(pid)
    t0 = time.time()
    ro = RunOptions(
        run_hydro=opts.run_hydro, run_flowlines=opts.run_flowlines,
        run_chemistry=opts.run_chemistry, run_reserves=opts.run_reserves,
        interpolation=opts.interpolation, variogram_model=opts.variogram_model,
        variogram_sill=opts.variogram_sill, variogram_range=opts.variogram_range,
        variogram_nugget=opts.variogram_nugget,
        anisotropy=tuple(opts.anisotropy) if len(opts.anisotropy) == 3 else (1.0, 1.0, 1.0),
        kv_ratio=opts.kv_ratio,
        stagnation_threshold=opts.stagnation_threshold,
        max_streamlines=opts.max_streamlines,
        contour_polygon=opts.contour_polygon,
        acid_threshold=opts.acid_threshold,
        output_every=opts.output_every,
        reactions=tuple(opts.reactions) if len(opts.reactions) == 4 else (True, True, True, True),
    )
    try:
        res = run_pipeline(p, ro)
    except Exception as e:  # noqa: BLE001
        import traceback
        raise HTTPException(500, f"Ошибка расчёта: {e}\n{traceback.format_exc()[-1500:]}")
    res["elapsed_sec"] = round(time.time() - t0, 2)
    return res


@app.get("/api/projects/{pid}/field")
def get_field(
    pid: str,
    name: str = Query(..., description="head | speed | k | uranium | porosity | acid | uranium_final | residual_final"),
    axis: str = Query("z", pattern="^(x|y|z)$"),
    index: int = Query(0),
    stride: int = Query(1, ge=1, le=10),
):
    p = _need(pid)
    arrays = p.arrays
    if name not in arrays:
        raise HTTPException(404, f"Поле '{name}' недоступно — выполните расчёт")
    arr = np.asarray(arrays[name], dtype=float)
    if axis == "z":
        if not (0 <= index < arr.shape[2]):
            raise HTTPException(400, "Индекс слоя вне диапазона")
        sl = arr[:, :, index]
        xaxis, yaxis = 0, 1
    elif axis == "y":
        sl = arr[:, index, :]
        xaxis, yaxis = 0, 2
    else:
        sl = arr[index, :, :]
        xaxis, yaxis = 1, 2
    sl = sl[::stride, ::stride]
    return {
        "name": name, "axis": axis, "index": index, "stride": stride,
        # плоский массив в порядке i * ny + j (i — индекс по X, j — по Y)
        "data": np.nan_to_num(sl, nan=0.0).ravel().tolist(),
        "min": float(np.nanmin(arr)), "max": float(np.nanmax(arr)),
        "shape": [int(sl.shape[0]), int(sl.shape[1])],
    }


@app.get("/api/projects/{pid}/profile")
def get_profile(pid: str, name: str = Query(...), i: int = Query(...), j: int = Query(...)):
    """Вертикальный профиль по колонке ячеек."""
    p = _need(pid)
    arrays = p.arrays
    if name not in arrays:
        raise HTTPException(404, "Поле недоступно")
    arr = np.asarray(arrays[name], dtype=float)
    if not (0 <= i < arr.shape[0] and 0 <= j < arr.shape[1]):
        raise HTTPException(400, "Индексы вне диапазона")
    return {
        "name": name, "i": i, "j": j,
        "z": [p.grid.z0 + (k + 0.5) * p.grid.dz for k in range(arr.shape[2])],
        "values": np.nan_to_num(arr[i, j, :]).tolist(),
    }


@app.get("/api/projects/{pid}/streamlines")
def get_streamlines(pid: str, limit: int = Query(300, ge=1, le=3000)):
    p = _need(pid)
    sls = p.results.get("streamlines") or []
    return {"count": len(sls), "streamlines": sls[:limit]}


@app.get("/api/projects/{pid}/flowline-chemistry")
def flowline_chemistry(pid: str, index: int = Query(0)):
    """Расчёт концентраций вдоль выбранной линии тока."""
    from .chemistry import simulate_flowline
    from .geo_ops import sample_field
    p = _need(pid)
    sls = p.results.get("_streamlines") or []
    if not sls:
        raise HTTPException(400, "Линии тока не рассчитаны")
    sl = sls[min(index, len(sls) - 1)]
    arrays = p.arrays
    grade = arrays.get("uranium")
    poro = arrays.get("porosity")
    pts = np.array(sl.points)
    if grade is None:
        raise HTTPException(400, "Нет модели урана")
    g_vals = np.array([sample_field(grade, p.grid, pt) for pt in pts])
    p_vals = np.array([sample_field(poro, p.grid, pt) if poro is not None else 0.15 for pt in pts])
    prof = simulate_flowline(pts, np.array(sl.times), g_vals, p_vals, p.chem)
    return {"index": index, "well_in": sl.well_in, "well_out": sl.well_out,
            "residence": sl.residence, **prof}


@app.get("/api/projects/{pid}/export")
def export_csv(pid: str, what: str = Query("wells")):
    p = _need(pid)
    if what == "wells":
        df = importers.wells_to_dataframe(p.wells)
    elif what == "properties":
        df = importers.properties_to_dataframe(p.props)
    elif what == "curves":
        c = p.results.get("curves")
        if not c:
            raise HTTPException(400, "Нет результатов расчёта")
        df = __import__("pandas").DataFrame(c)
    else:
        raise HTTPException(400, "Неизвестно что экспортировать")
    buf = io.StringIO()
    df.to_csv(buf, index=False, sep=";", encoding="utf-8-sig")
    return Response(
        content=buf.getvalue().encode("utf-8"),
        media_type="text/csv",
        headers={"Content-Disposition": f'attachment; filename="{what}.csv"'},
    )


@app.post("/api/projects/{pid}/save")
def save_project(pid: str):
    _need(pid)
    path = store.save(pid)
    return {"ok": True, "path": path}


@app.post("/api/projects/open")
async def open_project(file: UploadFile = File(...)):
    data = await file.read()
    path = os.path.join(store.base_dir, "_upload.s3dp")
    with open(path, "wb") as f:
        f.write(data)
    pid = store.load(path)
    os.remove(path)
    s = _project_summary(store.get(pid))
    s["id"] = pid
    return s


@app.api_route("/api/demo", methods=["GET", "POST"])
def load_demo():
    """Загрузка демонстрационного проекта (синтетический блок 5-2-1)."""
    from .demo import build_demo_project
    p = store.create("Демо: блок 5-2-1 (учебные данные)")
    build_demo_project(p)
    pid = _pid_of(p)
    s = _project_summary(p)
    s["id"] = pid
    return s


# ------------------------------------------------------------------ статика

if os.path.isdir(FRONTEND_DIR):
    app.mount("/", StaticFiles(directory=FRONTEND_DIR, html=True), name="frontend")
