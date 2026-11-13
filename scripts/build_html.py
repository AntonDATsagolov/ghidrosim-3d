#!/usr/bin/env python3
"""Сборка автономного HTML-файла «ГидроСимулятор 3D».

Результат — один файл, который открывается двойным щелчком в любом
браузере: Python, Node и сервер не нужны.

Два варианта сборки:
  * cdn     (по умолчанию) — three.js подключается с unpkg, размер файла
             около 230 КБ; требуется доступ в интернет;
  * offline — three.js встраивается в файл целиком (~2,2 МБ),
             интернет не требуется.

    python3 scripts/build_html.py                 # cdn
    python3 scripts/build_html.py --mode offline  # офлайн
"""

from __future__ import annotations

import argparse
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FRONTEND = os.path.join(ROOT, "frontend")
DIST = os.path.join(ROOT, "dist")

CDN_BASE = "https://unpkg.com/three@0.187.0"
IMPORTMAP_CDN = (
    '<script type="importmap">\n{\n  "imports": {\n'
    f'    "three": "{CDN_BASE}/build/three.module.js",\n'
    f'    "three/addons/": "{CDN_BASE}/examples/jsm/"\n'
    '  }\n}\n</script>\n'
)


def read(path: str) -> str:
    with open(path, encoding="utf-8") as f:
        return f.read()


# ------------------------------------------------------------------ обработчики кода

def strip_export(code: str) -> str:
    """Убирает завершающий export {...} — имена остаются в области видимости."""
    return re.sub(r"\nexport \{[^}]*\};\s*$", "\n", code)


def strip_exports(code: str) -> str:
    """Убирает export-объявления (в едином модуле они не нужны)."""
    code = re.sub(r"^\s*export\s+(?=(async\s+)?(function|class|const|let|var)\b)",
                  "", code, flags=re.MULTILINE)
    code = re.sub(r"^\s*export\s*\{[^}]*\};?\s*$", "", code, flags=re.MULTILINE)
    return code


def strip_import(code: str) -> str:
    """Убирает import-инструкции (всё встроено в один модуль)."""
    code = re.sub(r"^\s*import\s+[\s\S]*?from\s+'[^']+';\s*$", "", code, flags=re.MULTILINE)
    code = re.sub(r"^\s*import\s+'[^']+';\s*$", "", code, flags=re.MULTILINE)
    return code


def strip_reexport(code: str) -> str:
    """Убирает реэкспорт export { ... } from './three.core.js'."""
    return re.sub(r"export\s*\{[^}]*\}\s*from\s*'[^']+';", "", code)


def top_level_names(code: str) -> set[str]:
    """Имена, объявленные на верхнем уровне модуля."""
    names: set[str] = set()
    pattern = re.compile(
        r"^(?:const|let|var|function|class)\s+([A-Za-z_$][\w$]*)", re.MULTILINE)
    for m in pattern.finditer(code):
        names.add(m.group(1))
    for m in re.finditer(r"export\s*\{([^}]*)\}", code):
        for item in m.group(1).split(","):
            item = item.strip().split(" as ")[0].strip()
            if item and re.match(r"^[A-Za-z_$][\w$]*$", item):
                names.add(item)
    return names


def rename_conflicts(code: str, taken: set[str]) -> tuple[str, int]:
    """Переименовывает объявленные имена, уже занятые в three.core.js."""
    conflicts = sorted(top_level_names(code) & taken)
    for name in conflicts:
        new = name + "__mod"
        code = re.sub(r"\b" + re.escape(name) + r"\b", new, code)
        taken.add(new)
    return code, len(conflicts)


# ------------------------------------------------------------------ сборка

def app_code(mode: str) -> str:
    """Код приложения без внутренних импортов (всё в одном модуле).

    В режиме cdn импорты three.js сохраняются — они разрешаются через
    importmap; внутренние (относительные) импорты убираются всегда.
    """
    viewer = read(os.path.join(FRONTEND, "js", "viewer3d.js"))
    if mode == "cdn":
        viewer = viewer.replace("from '../vendor/OrbitControls.js'",
                                "from 'three/addons/controls/OrbitControls.js'")
    else:
        viewer = strip_import(viewer)
    parts = [
        viewer,
        read(os.path.join(FRONTEND, "js", "charts.js")),
        read(os.path.join(FRONTEND, "js", "engine.js")),
        read(os.path.join(FRONTEND, "js", "demo.js")),
        strip_import(read(os.path.join(FRONTEND, "js", "standalone.js"))),
    ]
    return "\n".join(parts)


def build(mode: str) -> str:
    html = read(os.path.join(FRONTEND, "index.html"))
    css = read(os.path.join(FRONTEND, "css", "style.css"))

    html = html.replace(
        '<link rel="stylesheet" href="css/style.css">',
        "<style>\n" + css + "\n</style>")
    html = html.replace('<button id="btnSave" class="btn btn-ghost">Сохранить</button>',
                        '<button id="btnSave" class="btn btn-ghost">Сохранить (JSON)</button>')
    html = html.replace('<button id="btnOpen" class="btn btn-ghost">Открыть</button>',
                        '<button id="btnOpen" class="btn btn-ghost">Открыть (JSON)</button>')
    html = html.replace(
        "<title>ГидроСимулятор 3D — гидродинамическое моделирование ПСВ урана</title>",
        "<title>ГидроСимулятор 3D — автономная версия (без Python)</title>")

    app = app_code(mode)

    # исходный importmap не нужен ни в одном из режимов
    html = re.sub(r'<script type="importmap">.*?</script>\s*', "", html, flags=re.DOTALL)

    if mode == "cdn":
        # three.js подключается с CDN: его имена живут в отдельном модуле,
        # поэтому собственные объявления приложения им не мешают.
        bundle = ("window.__ghidrosim3d = { version: '1.0', mode: 'cdn' };\n"
                  + strip_exports(app))
        html = html.replace('<script type="module" src="js/app.js"></script>',
                            IMPORTMAP_CDN + '<script type="module" src="js/app.js"></script>')
    else:
        # Офлайн: three.js встраивается целиком.
        three_core = strip_export(read(os.path.join(FRONTEND, "vendor", "three.core.js")))
        three_module = strip_exports(strip_import(
            strip_reexport(read(os.path.join(FRONTEND, "vendor", "three.module.js")))))
        three_module, _ = rename_conflicts(three_module, top_level_names(three_core))

        orbit = strip_import(read(os.path.join(FRONTEND, "vendor", "OrbitControls.js")))
        orbit = re.sub(r"\nexport \{[^}]*\};\s*$", "\n", orbit)
        orbit = "const OrbitControls = (function () {\n" + orbit + "\nreturn OrbitControls;\n})();"

        app = strip_exports(app).replace("THREE.", "")
        bundle = (
            "window.__ghidrosim3d = { version: '1.0', mode: 'offline' };\n"
            "/* === three.js r187 (MIT): ядро === */\n" + three_core +
            "\n/* === three.js r187 (MIT): WebGL-рендерер === */\n" + three_module +
            "\n/* === OrbitControls === */\n" + orbit +
            "\n/* === приложение === */\n{\n" + app + "\n}\n"
        )

    bundle = bundle.replace("</script", "<\\/script")

    html = re.sub(
        r'<script type="module" src="js/app\.js"></script>',
        lambda m: '<script type="module">\n' + bundle + "\n</" + "script>",
        html)

    os.makedirs(DIST, exist_ok=True)
    name = "ghidrosim3d.html" if mode == "cdn" else "ghidrosim3d-offline.html"
    out = os.path.join(DIST, name)
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Сборка автономного HTML")
    ap.add_argument("--mode", choices=["cdn", "offline"], default="cdn")
    args = ap.parse_args()
    out = build(args.mode)
    size = os.path.getsize(out)
    print(f"Собрано: {out}")
    print(f"Размер: {size/1024:.0f} КБ ({size:,} байт)")
    with open(out, encoding="utf-8") as f:
        content = f.read()
    for probe in ("ГидроСимулятор 3D", "buildPropertyModels", "ГидроСимулятор"):
        if probe not in content:
            print(f"ОШИБКА: в сборке нет фрагмента «{probe}»", file=sys.stderr)
            return 1
    print("Проверка пройдена.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
