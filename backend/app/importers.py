"""Импорт и экспорт данных скважин (CSV / XLSX) с автоматическим сопоставлением
столбцов и возможностью ручного уточнения.

Поддерживаемые поля:
    name       — наименование скважины
    x, y, z    — координаты устья (или отметка фильтра), м
    k          — коэффициент фильтрации, м/сут
    porosity   — пористость, д.е.
    uranium    — содержание урана, кг/т
    clay       — доля глинистой компоненты, 0..1
    lithology  — литология (текст)
    rate       — дебит, м3/сут
    well_type  — тип скважины (закачка/откачка)
    acid_in    — кислотность выщелачивающего раствора, г/л
    screen_top, screen_bottom — интервал фильтра, м
"""

from __future__ import annotations

import io
import re

import numpy as np
import pandas as pd

FIELD_ALIASES: dict[str, list[str]] = {
    "name": ["скважина", "наименование", "название", "well", "wellname", "well_name", "name",
             "№ скважины", "номер скважины", "скв", "скв.", "id", "borehole"],
    "x": ["x", "координата x", "координатах", "координата_х", "восток", "восточный", "coord_x",
          "coordinate_x", "easting"],
    "y": ["y", "координата y", "координатау", "координата_у", "север", "северный", "coord_y",
          "coordinate_y", "northing"],
    "z": ["z", "глубина", "отметка", "абс отметка", "абсолютная отметка", "отметка фильтра",
          "elev", "elevation", "coord_z", "глубина фильтра"],
    "k": ["кф", "коэффициент фильтрации", "коэффициент_фильтрации", "проницаемость", "k",
          "permeability", "kh", "кп", " filtration"],
    "porosity": ["пористость", "porosity", "poro", "м", "не пористость", "п"],
    "uranium": ["уран", "содержание урана", "cу", "cu", "u", "uranium", "grade", "запас",
                "содержание", "уран, кг/т", "u_pct"],
    "clay": ["глина", "доля глины", "глинистость", "clay", "глин", "песчаник"],
    "lithology": ["литология", "lithology", "lit", "тип породы", "порода", "lith"],
    "rate": ["дебит", "расход", "q", "rate", "производительность", "приток", "деб", "q_нагн"],
    "well_type": ["тип", "назначение", "type", "well_type", "категория", "скважина назначения"],
    "acid_in": ["кислотность", "кислота", "h2so4", "acid", "вр", "концентрация кислоты"],
    "screen_top": ["верхний интервал", "топ", "top", "начало интервала", "screen_top",
                   "интервал перфорации верх", "верх фильтра", "верхний фильтр",
                   "верх отверстия", "начало фильтра", "верх перфорации", "отметка верха фильтра"],
    "screen_bottom": ["нижний интервал", "низ", "bottom", "конец интервала", "screen_bottom",
                      "интервал перфорации низ", "низ фильтра", "нижний фильтр",
                      "низ отверстия", "конец фильтра", "низ перфорации", "отметка низа фильтра"],
}

NUMERIC_FIELDS = {"x", "y", "z", "k", "porosity", "uranium", "clay", "rate",
                  "acid_in", "screen_top", "screen_bottom"}


def _norm(s: str) -> str:
    return re.sub(r"[^a-zа-я0-9]+", "_", str(s).strip().lower()).strip("_")


def auto_map_columns(df: pd.DataFrame) -> dict[str, str | None]:
    """Сопоставление столбцов файла полям модели. Возвращает {поле: столбец}."""
    norm_cols = {_norm(c): c for c in df.columns}
    mapping: dict[str, str | None] = {}
    used: set[str] = set()
    for field, aliases in FIELD_ALIASES.items():
        found = None
        for alias in aliases:
            a = _norm(alias)
            if a in norm_cols and norm_cols[a] not in used:
                found = norm_cols[a]
                break
        if found is None:
            # нечёткое совпадение
            for nc, orig in norm_cols.items():
                if orig in used:
                    continue
                for alias in aliases:
                    a = _norm(alias)
                    if a and (a in nc or nc in a) and len(a) >= 3:
                        found = orig
                        break
                if found:
                    break
        if found:
            mapping[field] = found
            used.add(found)
    for field in FIELD_ALIASES:
        mapping.setdefault(field, None)
    return mapping


def read_table(data: bytes, filename: str) -> pd.DataFrame:
    """Чтение таблицы из CSV или XLSX/XLS."""
    name = filename.lower()
    if name.endswith((".xlsx", ".xlsm", ".xls")):
        return pd.read_excel(io.BytesIO(data), header=0, dtype=object)
    encodings = ["utf-8-sig", "utf-8", "cp1251", "koi8-r", "latin-1"]
    last_err = None
    for enc in encodings:
        try:
            return pd.read_csv(io.BytesIO(data), sep=None, engine="python",
                               encoding=enc, header=0, dtype=object)
        except Exception as e:  # noqa: BLE001
            last_err = e
    raise ValueError(f"Не удалось прочитать файл: {last_err}")


def _to_float(v) -> float:
    if v is None:
        return float("nan")
    if isinstance(v, (int, float, np.integer, np.floating)):
        return float(v)
    s = str(v).strip().replace(",", ".")
    s = re.sub(r"[^\d.+\-eE]", "", s)
    if s in ("", "-", "+", "."):
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def _norm_type(v) -> str:
    s = str(v).strip().lower()
    if s in ("закачка", "закачивающая", "нагн", "нагнетательная", "injection", "зак"):
        return "injection"
    if s in ("откачка", "откачивающая", "добыча", "отбор", "production", "отк"):
        return "production"
    if s in ("наблюдательная", "контрольная", "observation", "набл"):
        return "observation"
    return "injection"


def rows_to_wells(df: pd.DataFrame, mapping: dict[str, str | None]) -> list[dict]:
    """Преобразование таблицы в список записей скважин."""
    out: list[dict] = []
    for _, row in df.iterrows():
        rec: dict = {}
        for field, col in mapping.items():
            if not col or col not in df.columns:
                rec[field] = None
                continue
            val = row[col]
            if field in NUMERIC_FIELDS:
                rec[field] = _to_float(val)
            elif field == "lithology":
                rec[field] = None if pd.isna(val) else str(val).strip()
            elif field == "well_type":
                rec[field] = _norm_type(val)
            else:
                rec[field] = None if pd.isna(val) else str(val).strip()
        if rec.get("x") is None or rec.get("y") is None:
            continue
        out.append(rec)
    return out


def wells_to_dataframe(wells) -> pd.DataFrame:
    """Экспорт скважин проекта в таблицу."""
    rows = []
    for w in wells:
        x, y = w.head_xy()
        z = w.screen.top if w.screen else 0.0
        rows.append({
            "name": w.name,
            "x": x,
            "y": y,
            "z": z,
            "well_type": w.well_type.value,
            "rate": w.rate,
            "acid_in": w.acid_in,
            "screen_top": w.screen.top if w.screen else None,
            "screen_bottom": w.screen.bottom if w.screen else None,
            "notes": w.notes,
        })
    return pd.DataFrame(rows)


def properties_to_dataframe(props) -> pd.DataFrame:
    """Экспорт скважинных данных свойств в таблицу."""
    rows = []
    for x, y, z, v in props.permeability:
        rows.append({"x": x, "y": y, "z": z, "k": v})
    for x, y, z, v in props.uranium:
        rows.append({"x": x, "y": y, "z": z, "uranium": v})
    for x, y, z, v in props.porosity:
        rows.append({"x": x, "y": y, "z": z, "porosity": v})
    for x, y, z, v in props.clay:
        rows.append({"x": x, "y": y, "z": z, "clay": v})
    for x, y, z, s in props.lithology:
        rows.append({"x": x, "y": y, "z": z, "lithology": s})
    return pd.DataFrame(rows)
