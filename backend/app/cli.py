"""Точка входа: запуск сервера «ГидроСимулятор 3D».

    python -m app.cli            # http://127.0.0.1:8000
    python -m app.cli --port 9000 --reload
"""

from __future__ import annotations

import argparse
import os
import sys
import webbrowser


def main() -> int:
    ap = argparse.ArgumentParser(description="ГидроСимулятор 3D")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--reload", action="store_true")
    ap.add_argument("--no-browser", action="store_true")
    args = ap.parse_args()

    app_dir = os.path.dirname(os.path.abspath(__file__))        # .../backend/app
    backend_dir = os.path.dirname(app_dir)                       # .../backend
    root = os.path.dirname(backend_dir)                          # корень проекта
    for p in (backend_dir, root):
        if p not in sys.path:
            sys.path.insert(0, p)
    os.chdir(root)

    try:
        import uvicorn
    except ImportError:
        print("Не установлен uvicorn: pip install -r requirements.txt", file=sys.stderr)
        return 1

    url = f"http://{args.host}:{args.port}"
    print(f"ГидроСимулятор 3D → {url}")
    if not args.no_browser:
        try:
            webbrowser.open(url)
        except Exception:  # noqa: BLE001
            pass

    uvicorn.run(
        "app.main:app",
        host=args.host,
        port=args.port,
        reload=args.reload,
        log_level="info",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
