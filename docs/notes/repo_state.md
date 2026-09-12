# Состояние репозитория — что готово и что делать дальше

Обновляется в конце каждой сессии. Быстрый вход для новой сессии:
прочитать AGENTS.md (правила), затем этот файл, затем PLAN §8.

## Готово

Сессия 1 (12.09.2026): ядро oracle-1.0 (`oracle/`), BFS (`exact/`),
сертифицированные таблицы 4 <= n <= 12 (n=11, 12 не в git — воспроизведение
описано в `data/tables/dist_n1{1,2}.json`), C2v–C14v VERIFIED,
`bounds/known.py`, отчёт `data/runs/session1_report.md`.

Сессия 2 (12.09.2026): входные файлы добавлены владельцем в `docs/incoming/`
(N1 + прежний комплект; `verify_lrx.py`, `verification.json`,
`verify_lrx_gap_bound.py` отсутствуют).
- Аудит N1 завершён: `docs/notes/strict_proof_audit.md`; C16–C21, C23 PROVED
  (локальный аудит); C22 CONJECTURED без контрпримеров.
- Реализация конструкции по тексту: `constructions/strict_upper.py`;
  проверки `checks/check_C16.py`, отчёт `data/runs/strict_audit/report.md`.
- C9 доказано явным словом: `docs/proofs/C9_sigma_upper.md`, `checks/check_C9.py`.
- Удостоверенная верхняя граница при всех n >= 4: U_n (C16), `bounds/known.py`.

## Чего нет

- Аудит текста C1 и gap-файла (C4/C10/C11) — CLAIMED; `checks/check_C1.py`,
  `check_gap_T1.py` не написаны (PLAN §9.9). Скрипт C1 автора только запущен.
- Аудит S1 (нижняя оценка C8) — не начат; `docs/notes/literature_audit.md` нет.
- Общая верхняя оценка вида B_n + o(n) — нет; C15 OPEN.
- Модули §9.7 (invariants), §9.8 (evaluation) — не созданы.

## Следующая задача

PLAN §8 / §4.1: полная карта потерь конструкции N1 на худших входах n = 7, 8
и выборках n = 20, 50, 100: разложить `min_c[2F_c - S_c + H(c)] - d(pi)`
на (минимум против среднего по c), (H(c) против маршрута по переносчикам),
(сокращения на стыках); выбрать одно направление §4.2. Бюджет — 1 сессия.
Инструменты готовы: `strict_upper.shift_word`/`perm_summary` дают все
статистики правила 15; таблицы d(pi) — до n = 12.
