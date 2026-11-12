#!/bin/bash
# Запуск «ГидроСимулятор 3D»
cd "$(dirname "$0")" || exit 1

if [ ! -x ".venv/bin/python" ]; then
  echo "Создаём виртуальное окружение (первый запуск)…"
  python3 -m venv .venv || { echo "Не удалось создать .venv"; read -r -p "Enter…"; exit 1; }
  .venv/bin/pip install -q --upgrade pip
  .venv/bin/pip install -q -r requirements.txt || {
    echo "Не удалось установить зависимости"; read -r -p "Enter…"; exit 1; }
fi

# Сторонняя библиотека three.js (3D-вид) скачивается при первом запуске
if [ ! -f "frontend/vendor/three.module.js" ]; then
  echo "Загружаем three.js…"
  python3 scripts/fetch_vendor.py || echo "three.js не загружена — 3D-вид будет недоступен"
fi

echo "Запуск сервера: http://127.0.0.1:8000"
.venv/bin/python backend/app/cli.py --port 8000
