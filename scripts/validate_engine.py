#!/usr/bin/env python3
"""Сверка порта на JavaScript с эталонной Python-реализацией.

    python3 scripts/validate_engine.py

Скрипт выгружает демо-проект в JSON, считает его эталонным способом
(Python), запускает Node-версию на тех же данных и сравнивает результаты.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "backend"))

import numpy as np  # noqa: E402

from app import hydro  # noqa: E402
from app.demo import build_demo_project  # noqa: E402
from app.models import BoundaryCondition, Project  # noqa: E402
from app.project import RunOptions, run_pipeline  # noqa: E402

OPTIONS = {
    "interpolation": "kriging",
    "variogram_model": "spherical",
    "variogram_sill": 1.0,
    "variogram_range": 150.0,
    "variogram_nugget": 0.0,
    "kv_ratio": 0.2,
    "output_every": 30.0,
    "acid_threshold": 5.0,
    "stagnation_threshold": 0.05,
}


def dump_project(path: str) -> None:
    p = Project(name="validate")
    build_demo_project(p)
    faces = {f: (p.bc_faces.get(f) or "neumann") for f in hydro.FACE_AXES}
    data = {
        "project": {
            "grid": {
                "x0": p.grid.x0, "y0": p.grid.y0, "z0": p.grid.z0,
                "nx": p.grid.nx, "ny": p.grid.ny, "nz": p.grid.nz,
                "dx": p.grid.dx, "dy": p.grid.dy, "dz": p.grid.dz,
            },
            "wells": [
                {
                    "id": w.id, "name": w.name, "well_type": w.well_type.value,
                    "trajectory": [{"x": t.x, "y": t.y, "z": t.z} for t in w.trajectory],
                    "screens": [{"top": s.top, "bottom": s.bottom} for s in w.screens],
                    "rate": w.rate, "acid_in": w.acid_in,
                }
                for w in p.wells
            ],
            "props": {
                "permeability": p.props.permeability,
                "porosity": p.props.porosity,
                "uranium": p.props.uranium,
                "clay": p.props.clay,
                "lithology": p.props.lithology,
            },
            "chem": p.chem.__dict__,
            "bcFaces": faces,
            "bcHeads": p.bc_heads,
        },
        "options": OPTIONS,
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False)


def python_reference() -> dict:
    p = Project(name="validate")
    build_demo_project(p)
    res = run_pipeline(p, RunOptions(
        interpolation=OPTIONS["interpolation"],
        variogram_model=OPTIONS["variogram_model"],
        variogram_sill=OPTIONS["variogram_sill"],
        variogram_range=OPTIONS["variogram_range"],
        variogram_nugget=OPTIONS["variogram_nugget"],
        kv_ratio=OPTIONS["kv_ratio"],
        output_every=OPTIONS["output_every"],
        acid_threshold=OPTIONS["acid_threshold"],
        stagnation_threshold=OPTIONS["stagnation_threshold"],
    ))
    c = res["curves"]
    r = res["reserves"]
    return {
        "hydro": {
            "residual": res["hydro"]["residual"],
            "total_in": res["hydro"]["total_in"],
            "total_out": res["hydro"]["total_out"],
            "head_min": res["hydro"]["head_min"],
            "head_max": res["hydro"]["head_max"],
            "speed_max": res["hydro"]["speed_max"],
        },
        "property_stats": res["property_stats"],
        "streamlines": {
            "count": len(res.get("streamlines") or []),
            "reached": sum(1 for s in (res.get("streamlines") or []) if s["well_out"]),
        },
        "chemistry": {
            "times": c["time"],
            "pr_uranium": c["pr_uranium"],
            "vr_acid": c["vr_acid"],
            "cumulative_u": c["cumulative_u"],
            "cumulative_acid": c["cumulative_acid"],
            "prod_volume": c["prod_volume"],
        },
        "reserves": {
            "reserve_kg": r["reserve_kg"],
            "reserve_t": r["reserve_t"],
            "acidified_cells": r["acidified_cells"],
            "acidified_volume": r["acidified_volume"],
            "area_ha": r["area_ha"],
            "avg_grade": r["avg_grade"],
            "avg_thickness": r["avg_thickness"],
            "weighted_metre_percent": r["weighted_metre_percent"],
        },
    }


def compare(name: str, a, b, rel: float = 2e-3, abs_: float = 1e-6) -> bool:
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        scale = max(abs(a), abs(b), abs_)
        diff = abs(a - b) / scale
        ok = diff <= rel
        print(f"  {'OK ' if ok else 'РАСХОЖДЕНИЕ'} {name}: python={a:.6g} js={b:.6g} "
              f"(отн. откл. {diff:.2e})")
        return ok
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            print(f"  РАСХОЖДЕНИЕ {name}: разная длина {len(a)} vs {len(b)}")
            return False
        worst = 0.0
        worst_i = -1
        for i, (x, y) in enumerate(zip(a, b)):
            scale = max(abs(x), abs(y), abs_)
            d = abs(x - y) / scale
            if d > worst:
                worst, worst_i = d, i
        ok = worst <= rel
        print(f"  {'OK ' if ok else 'РАСХОЖДЕНИЕ'} {name}: {len(a)} значений, "
              f"макс. откл. {worst:.2e}" + (f" (индекс {worst_i})" if worst_i >= 0 else ""))
        return ok
    ok = a == b
    print(f"  {'OK ' if ok else 'РАСХОЖДЕНИЕ'} {name}: python={a} js={b}")
    return ok


def main() -> int:
    node = "node"
    with tempfile.TemporaryDirectory() as tmp:
        inp = os.path.join(tmp, "input.json")
        outp = os.path.join(tmp, "js.json")
        dump_project(inp)
        proc = subprocess.run([node, os.path.join(ROOT, "scripts", "validate_engine.mjs"), inp, outp],
                              capture_output=True, text=True, timeout=900)
        if proc.returncode != 0:
            print("Ошибка запуска Node:", proc.stdout, proc.stderr)
            return 1
        print(proc.stdout.strip())
        with open(outp, encoding="utf-8") as f:
            js = json.load(f)

    py = python_reference()
    print("\nСравнение Python-эталона и JS-порта:")
    ok = True
    # невязка оценивается относительно суммарного расхода (абсолютное
    # значение в м3/сут слишком мало для относительного критерия)
    py_rel = abs(py["hydro"]["residual"]) / max(py["hydro"]["total_in"], 1.0)
    js_rel = abs(js["hydro"]["residual"]) / max(js["hydro"]["total_in"], 1.0)
    bal_ok = py_rel < 1e-9 and js_rel < 1e-9
    print(f"  {'OK ' if bal_ok else 'РАСХОЖДЕНИЕ'} гидродинамика: невязка (отн. к расходу): "
          f"python={py_rel:.2e} js={js_rel:.2e} (порог 1e-09)")
    ok &= bal_ok
    ok &= compare("гидродинамика: закачка", py["hydro"]["total_in"], js["hydro"]["total_in"])
    ok &= compare("гидродинамика: отбор", py["hydro"]["total_out"], js["hydro"]["total_out"])
    ok &= compare("гидродинамика: напор min", py["hydro"]["head_min"], js["hydro"]["head_min"])
    ok &= compare("гидродинамика: напор max", py["hydro"]["head_max"], js["hydro"]["head_max"])
    ok &= compare("гидродинамика: скорость max", py["hydro"]["speed_max"], js["hydro"]["speed_max"])
    ok &= compare("свойства: k min", py["property_stats"]["k_min"], js["property"]["k_min"])
    ok &= compare("свойства: k max", py["property_stats"]["k_max"], js["property"]["k_max"])
    ok &= compare("свойства: k среднее", py["property_stats"]["k_mean"], js["property"]["k_mean"])
    ok &= compare("свойства: U среднее", py["property_stats"]["u_mean"], js["property"]["u_mean"])
    ok &= compare("линии тока: количество", py["streamlines"]["count"], js["streamlines"]["count"], rel=0.05)
    ok &= compare("линии тока: дошли до отбора", py["streamlines"]["reached"], js["streamlines"]["reached"], rel=0.1)
    ok &= compare("ПР (мг/л)", py["chemistry"]["pr_uranium"], js["chemistry"]["pr_uranium"])
    ok &= compare("ВР (г/л)", py["chemistry"]["vr_acid"], js["chemistry"]["vr_acid"])
    ok &= compare("накоплено U (кг)", py["chemistry"]["cumulative_u"], js["chemistry"]["cumulative_u"])
    ok &= compare("расход кислоты (кг)", py["chemistry"]["cumulative_acid"], js["chemistry"]["cumulative_acid"])
    ok &= compare("объём отбора (м3)", py["chemistry"]["prod_volume"], js["chemistry"]["prod_volume"])
    ok &= compare("запасы (т)", py["reserves"]["reserve_t"], js["reserves"]["reserve_t"], rel=5e-3)
    ok &= compare("закисляемых ячеек", py["reserves"]["acidified_cells"], js["reserves"]["acidified_cells"], rel=5e-3)
    ok &= compare("ср. содержание (кг/т)", py["reserves"]["avg_grade"], js["reserves"]["avg_grade"])
    ok &= compare("ср. мощность (м)", py["reserves"]["avg_thickness"], js["reserves"]["avg_thickness"])
    ok &= compare("метропроцент", py["reserves"]["weighted_metre_percent"], js["reserves"]["weighted_metre_percent"])

    print("\nИтог:", "реализации совпадают" if ok else "ЕСТЬ РАСХОЖДЕНИЯ")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
