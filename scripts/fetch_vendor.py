#!/usr/bin/env python3
"""Загрузка сторонней библиотеки three.js в frontend/vendor.

three.js (r187) распространяется под лицензией MIT и не входит в состав
репозитория: файлы проверяются по SHA-256 и скачиваются при первом запуске.
Так репозиторий остаётся лёгким, а программа работает и без интернета
после однократной загрузки.

    python3 scripts/fetch_vendor.py [--force]
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
VENDOR = ROOT / "frontend" / "vendor"
BASE = "https://unpkg.com/three@0.187.0"

FILES = [
    {
        "path": "three.module.js",
        "url": f"{BASE}/build/three.module.js",
        "sha256": "104f56145899722094556490d787aae0617cc2fc3382f31df02d724ea1aee341",
    },
    {
        "path": "three.core.js",
        "url": f"{BASE}/build/three.core.js",
        "sha256": "613e370f76638fdef83926642bf256dc9615684922de1ce1f63a8583e5d78578",
    },
    {
        "path": "OrbitControls.js",
        "url": f"{BASE}/examples/jsm/controls/OrbitControls.js",
        "sha256": "e7fca2fe9e7670a828af9b28084669a5909445b7f36f1f26cee9bf8f2961c32c",
    },
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(item: dict, force: bool = False) -> bool:
    target = VENDOR / item["path"]
    if target.exists() and not force:
        if sha256(target) == item["sha256"]:
            print(f"[ok]   {item['path']} уже загружен")
            return True
        print(f"[warn] {item['path']}: контрольная сумма не совпала, перезагрузка")
    print(f"[get]  {item['url']}")
    with urllib.request.urlopen(item["url"], timeout=60) as resp:
        data = resp.read()
    digest = hashlib.sha256(data).hexdigest()
    if digest != item["sha256"]:
        print(f"[fail] {item['path']}: SHA-256 не совпал\n"
              f"       ожидалось {item['sha256']}\n"
              f"       получено  {digest}", file=sys.stderr)
        return False
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)
    print(f"[done] {item['path']} ({len(data):,} байт)")
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description="Загрузка three.js в frontend/vendor")
    ap.add_argument("--force", action="store_true", help="перезагрузить даже если файлы есть")
    args = ap.parse_args()
    ok = all(fetch(item, args.force) for item in FILES)
    if not ok:
        print("Не удалось загрузить three.js — 3D-вид будет недоступен.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
