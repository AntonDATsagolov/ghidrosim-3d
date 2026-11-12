"""Тесты геостатистики, запасов, импорта и полного конвейера."""

from __future__ import annotations

import numpy as np
import pytest

from app import geostat, reserves
from app.demo import build_demo_project
from app.models import GridSpec, Project
from app.project import RunOptions, run_pipeline


# ---------------------------------------------------------------- геостатистика

def test_idw_exact_at_data_point():
    pts = np.array([[0.0, 0.0, 0.0], [100.0, 0.0, 0.0], [0.0, 100.0, 0.0]])
    vals = np.array([1.0, 2.0, 3.0])
    s = geostat.GeoSettings(method="idw")
    out = geostat.interpolate(pts, vals, pts, s)
    assert np.allclose(out, vals, atol=1e-6)


def test_kriging_interpolates_data_points():
    rng = np.random.default_rng(0)
    pts = rng.random((25, 3)) * 100
    vals = pts[:, 0] * 2 + pts[:, 1]
    s = geostat.GeoSettings(method="kriging",
                            variogram=geostat.Variogram(model="spherical", sill=1.0, range_=80.0))
    out = geostat.interpolate(pts, vals, pts, s)
    # кригинг должен воспроизводить данные с точностью порядка анизотропии ряда
    assert np.corrcoef(out, vals)[0, 1] > 0.97
    assert np.abs(out - vals).max() < 0.25 * vals.std()


def test_variogram_models_are_monotonic():
    for model in ("spherical", "exponential", "gaussian", "linear"):
        vg = geostat.Variogram(model=model, sill=2.0, range_=50.0)
        h = np.linspace(0, 120, 50)
        g = vg.evaluate(h)
        assert np.all(np.diff(g) >= -1e-9)
        assert g[0] == pytest.approx(0.0)
        assert vg.evaluate(1e6) == pytest.approx(2.0, abs=0.6)


def test_nearest_lithology():
    pts = np.array([[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]])
    labels = ["песчаник", "алевролит"]
    out = geostat.nearest_lithology(pts, labels, np.array([[5.0, 0.0, 0.0], [95.0, 0.0, 0.0]]))
    assert list(out) == ["песчаник", "алевролит"]


# ---------------------------------------------------------------- запасы

def test_reserve_calculation():
    grid = GridSpec(nx=4, ny=4, nz=2, dx=10.0, dy=10.0, dz=5.0)
    grade = np.full(grid.shape, 1.0)
    acid = np.zeros(grid.shape)
    acid[:, :, 0] = 20.0     # закислен только верхний слой
    acid[:, :, 1] = 1.0
    res = reserves.calculate(grid, grade, acid, acid_threshold=5.0)
    assert res.acidified_cells == 16
    assert res.acidified_volume == 16 * grid.cell_volume
    # запас: 16 ячеек по 1 кг/т * масса руды в ячейке
    cell_ore_t = 2200 * grid.cell_volume / 1000.0
    assert res.reserve_kg == pytest.approx(16 * 1.0 * cell_ore_t, rel=1e-6)
    assert res.avg_grade == pytest.approx(1.0)
    assert res.avg_thickness == pytest.approx(5.0)


def test_contour_polygon():
    grid = GridSpec(nx=10, ny=10, nz=2, dx=10.0, dy=10.0, dz=5.0)
    poly = [[20.0, 20.0], [80.0, 20.0], [80.0, 80.0], [20.0, 80.0]]
    mask = reserves.contour_from_polygon(grid, poly)
    assert mask.any()
    assert not mask.all()
    # все выбранные ячейки внутри полигона
    xc, yc, _ = grid.centers()
    assert np.all((xc[mask] > 15) & (xc[mask] < 85))
    assert np.all((yc[mask] > 15) & (yc[mask] < 85))


# ---------------------------------------------------------------- полный конвейер

def test_demo_pipeline_runs_and_is_balanced():
    p = Project(name="test")
    build_demo_project(p)
    res = run_pipeline(p, RunOptions(output_every=90.0))
    assert res["hydro"]["balance_ok"]
    assert abs(res["hydro"]["residual"]) < 1e-3
    assert res["property_stats"]["k_min"] > 0
    assert res["reserves"]["reserve_kg"] > 0
    assert res["reserves"]["acidified_cells"] > 0
    curves = res["curves"]
    assert len(curves["time"]) == len(curves["pr_uranium"])
    assert curves["cumulative_u"][-1] > 0
    # концентрация урана в ПР — физически разумный диапазон
    assert 1.0 < max(curves["pr_uranium"]) < 5000.0
    # остаточный запас урана в руде снижается
    assert res["reserves"]["avg_grade"] <= res["property_stats"]["u_mean"] * 1.2


def test_demo_has_balanced_wells_inside_block():
    p = Project(name="test")
    build_demo_project(p)
    g = p.grid
    for w in p.wells:
        x, y = w.head_xy()
        assert g.x0 <= x <= g.x0 + g.nx * g.dx
        assert g.y0 <= y <= g.y0 + g.ny * g.dy
        assert w.screen_cells(g), f"скважина {w.name} вне сетки"
    inj = sum(w.rate for w in p.wells if w.well_type.value == "injection")
    prod = sum(w.rate for w in p.wells if w.well_type.value == "production")
    assert inj + prod == pytest.approx(0.0, abs=1e-6)


# ---------------------------------------------------------------- импорт

def test_import_auto_mapping_and_rows(tmp_path):
    from app import importers
    csv = (
        "Скважина;X;Y;Коэффициент фильтрации;Уран;Дебит\n"
        "1;100;100;2,5;0,8;500\n"
        "2;200;100;3,1;0,9;-480\n"
    )
    path = tmp_path / "wells.csv"
    path.write_text(csv, encoding="utf-8")
    df = importers.read_table(path.read_bytes(), "wells.csv")
    mapping = importers.auto_map_columns(df)
    assert mapping["x"] == "X"
    assert mapping["k"] == "Коэффициент фильтрации"
    rows = importers.rows_to_wells(df, mapping)
    assert len(rows) == 2
    assert rows[0]["x"] == 100.0
    assert rows[0]["k"] == pytest.approx(2.5)   # запятая как десятичный разделитель
    assert rows[1]["rate"] == pytest.approx(-480.0)
