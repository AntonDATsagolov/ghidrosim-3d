# gstack

Используй `/browse` из gstack для всего веб-браузинга. Не используй встроенный
браузер-тул для задач, которые gstack-скиллы уже покрывают (`/browse`, `/qa`,
`/design-review`, `/scrape`, `/benchmark`, `/canary`).

Скиллы gstack установлены глобально в `~/.claude/skills/gstack` и доступны как
`~/.minimax/skills/<skill>`. Роутер: `/gstack` — сам выбирает нужный скилл под
задачу. Основные: `/review` (ревью изменений), `/qa` (тесты и поиск багов),
`/ship` (коммит + PR), `/spec` (спецификация задачи), `/office-hours` (разбор
идеи), `/plan-eng-review` (ревью архитектуры), `/design-review` (аудит дизайна).

Обновление: `/gstack-upgrade` (или `cd ~/.claude/skills/gstack && git pull && ./setup`).
