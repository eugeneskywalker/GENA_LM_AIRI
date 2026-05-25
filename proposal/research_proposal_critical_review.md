# Critical Review of Research Proposal

**Объект ревью:** *GENA-LM: A Critical Review with Layer-Aware Reproducibility Experiments* (AIRI Summer School 2026 Research Proposal, автор: Korshunov E.)
**Дата ревью:** 2026-05-25
**Ревьюер:** строгий академический методолог
**Формат исходного proposal:** не типичный «новый research proposal», а **paper-critique-and-extension format** (заявка в школу AIRI). Это влияет на применимость стандартных критериев — некоторые блоки (например, «Research Questions», «Limitations») в исходнике отсутствуют, что и составляет часть критики.

---

## 1. Executive Summary

**Общее впечатление.** Proposal демонстрирует **сильную методологическую компетенцию для уровня bachelor-application**: автор воспроизвёл и провёл пять frozen-model экспериментов на V100, получил численные результаты с figures и таблицами, и связал их с актуальным литературным контекстом (DART-Eval, BERTology). Это редкое для уровня application качество исполнения.

**Главная сила:** реальный, воспроизводимый эмпирический материал. Не пересказ статьи, а **собственные эксперименты** с конкретными числами и графиками — большая редкость в AIRI заявках.

**Главная слабость:** **систематическое overstating** результатов. Несколько центральных claims (особенно «direct quantitative confirmation of DART-Eval», «causal evidence», «mid-layer geometric primacy») сформулированы строже, чем позволяют данные. Опытный ревьюер AIRI это заметит и downgrade оценку.

**Вторая по серьёзности слабость:** **отсутствие явного указания на число random seeds** для самых важных экспериментов (E3, E6). Главный сравнительный claim («HyenaDNA побеждает GENA-LM 3 из 4») висит на разницах в F1/AUROC порядка 0.01–0.05, что может быть в пределах seed noise.

**Можно ли подавать в текущем виде:** да — но **4–5 часов работы** (см. раздел 6) дадут максимальный прирост убедительности и снимут главные риски critique.

**Срочные исправления первого приоритета:**
1. Перезапустить E3 + E6 с n=3 seeds, добавить std в таблицы.
2. Снять формулировку «DART-Eval» → «DART-Eval-inspired».
3. Добавить qualifier «frozen-probing» ко всем claims про HyenaDNA winning.
4. Добавить минимум один negative control в E5 (мутации вне мотива).
5. Уточнить какая именно версия HyenaDNA и какое context length использовались.
6. Добавить отсутствующий блок **Limitations**.

---

## 2. Overall Assessment

| Критерий | Оценка | Комментарий |
|---|---|---|
| Ясность research idea | **Strong** | Понятно за 30 секунд: «layer-aware critique GENA-LM + head-to-head с альтернативой» |
| Логичность структуры | **Strong** | Strengths → Weaknesses → Experiments → Results — каноническая структура для critique proposal, аккуратно реализована |
| Убедительность research gap | **Acceptable but needs revision** | Gap определён через 5 weaknesses статьи (W3, W6, W7, W10, W12) — конкретно. Но gap «no DART-Eval evaluation» подменяется суррогатом, который не есть DART-Eval |
| Соответствие методологии целям | **Acceptable but needs revision** | Frozen probing — оправданный выбор для масштаба 1 студент × V100 × 2 недели. Но выводы делаются как для fine-tuned setting |
| Статистическая строгость | **Weak** | Seeds указаны только для E1; нет CI/std для большинства центральных чисел; нет multiple comparison correction; нет statistical tests |
| Готовность к подаче | **Acceptable but needs revision** | Подавать можно, но 5 часов работы превратят это из «middle of pack» в «strong submission» |

**Итоговый balanced verdict:** **Acceptable but needs revision.** Базовая работа сильная, но overstating + отсутствие seeds для главных claims — это прямые мишени критики, которые легко закрываются.

---

## 3. Major Issues

| # | Problem | Why It Matters | Severity | How to Fix |
|---|---|---|---|---|
| 1 | E6 называется «DART-Eval-style 4-way benchmark», но реально запущена задача из Genomic Benchmarks (`human_ocr_ensembl`), а не DART-Eval Task 1–5 | Главный «новый» вклад заявлен как DART-Eval validation. Если ревьюер знает DART-Eval, он увидит подмену и потеряет доверие к остальным claims | **Critical** | Переименовать в «DART-Eval-inspired calibration on a regulatory-DNA proxy task». В тексте: «we replicate the spirit of the DART-Eval critique on an OCR classification proxy; the formal DART-Eval pipeline is in follow-ups» |
| 2 | Для E3 и E6 не указано число seeds (для E1 явно сказано «n=3 seeds») | E3 и E6 — основа главного claim «HyenaDNA побеждает GENA-LM, подтверждая DART-Eval». Разницы 0.013 F1 на promoter и 0.030 AUROC на OCR могут быть в пределах seed noise. Без n≥3 + std весь central claim висит в воздухе | **Critical** | Перезапустить E3 и E6 с 3 seeds, добавить столбец std или CI в таблицы |
| 3 | Сравнение HyenaDNA vs GENA-LM делается на **frozen probing**, но интерпретируется как сравнение моделей в целом — авторские числа GENA-LM из NAR 2025 получены на **fine-tuned** setup | Сравнение режимов «как fine-tuned» vs «как frozen» — apples to oranges. Любой ревьюер из AIRI Bioinformatics (а это именно тот ревьюер, который и читает) знает что fine-tuned GENA-LM на promoter даёт 94.16% F1, а не 0.777 | **Critical** | Везде где есть claim про winning HyenaDNA — qualifier «in the frozen-probing regime». Дополнительно: одна fine-tuned точка (промоутер, 30 мин на V100) колоссально усиливает proposal |
| 4 | E5 «causal» claim r=0.893 без negative controls | Высокий r может объясняться tokenization artefacts BPE при OOD-нуклеотидах, а не «GENA-LM учит CTCF grammar». Без контроля (мутации вне мотива → r ≈ 0) claim не отделим от artefact-объяснения | **High** | Добавить минимум 200 forward passes: мутации в случайных позициях flank вне мотива. Ожидаемый r ≈ 0 (если r=0.893 валиден) или r ≈ 0.5+ (если artefact) |
| 5 | E1 «mid-layer = biology» — справедливо только для 2 из 4 задач | Coding peak L5 ✅, Promoter peak L4 ✅, **Enhancer peak L1**, **OCR peak L2**. Последние две в области chance (F1 ≈ 0.6–0.66), их «пик» — это шум. Generalisation «peak in middle» — overreach | **High** | Переформулировать: «for tasks the model meaningfully resolves (promoter, coding), peak is mid-layer; long-range regulatory tasks (enhancer, OCR) are near-chance regardless of layer, suggesting BigBird/RMT variants are required» |
| 6 | E4 cherry-picking метрики: Silhouette → L4, NMI → L0, ARI → L12 | Три метрики — три «лучших» слоя. Вывод «mid-layer best for unsupervised» делается по одной метрике (silhouette), при этом две другие говорят противоположное | **High** | Честно перечислить все три и формулировать: «mid-layer is best for geometric arrangement (silhouette); total information (NMI) and categorical structure (ARI) do not single out mid-layer» |
| 7 | E3 берёт last-layer GENA-LM (L12), хотя E1 показал что best layer = L4–L5 | Это buy ourselves vulnerability: ревьюер скажет «вы намеренно использовали худший слой для центрального claim». Если в E3 переключить на L4 — может быть GENA-LM не уступит | **High** | Перезапустить E3 с GENA-LM на best-per-task layer (из E1). Показать обе версии. Это закрывает риск и даёт дополнительный sound-bite («even at best layer, GENA-LM matches but does not beat HyenaDNA») |
| 8 | Limitations block отсутствует полностью | Любой proposal должен иметь Limitations. Отсутствие читается как «автор не понимает границ своей работы» | **High** | Добавить раздел Limitations (см. секцию 4.12 этого ревью) |
| 9 | DNABERT-2 «failed to load due to einops/transformers interaction» | Это технический сбой, не научный аргумент. Для ревьюера читается как «не разобрался с зависимостями». Слабит доверие к остальной технической работе | **Medium** | Либо починить (pin transformers версии), либо убрать упоминание полностью. Не оставлять в текущем виде |
| 10 | Не указано какая именно версия HyenaDNA и context length тестирования | `hyenadna-tiny-1k-seqlen` vs `hyenadna-tiny-16k-seqlen` — разные модели. Без явного указания ревьюер заподозрит cherry-picking | **Medium** | Добавить в таблицу setup: точный HF ID модели + max sequence length для обеих моделей |
| 11 | Mean-pooling BPE (≈110 tokens, медиана 9 bp/token) vs single-nuc (1000 positions) — несимметричная операция | Сравнение HyenaDNA vs GENA-LM может частично объясняться pooling strategy, не моделями. Sensitivity к pooling не проверена | **Medium** | Минимум: для одной задачи показать CLS + mean + max pooling, демонстрируя что вывод стабилен |
| 12 | Promoter «nontata» subset — почему именно non-tata? | Authors GENA-LM benchmarks on EPDnew (с TATA). Non-tata — более сложный subset → искусственно занижает F1. Без обоснования выглядит cherry-picked в обратную сторону | **Medium** | Либо обосновать выбор non-tata («harder subset, tests long-range context»), либо запустить на полном EPDnew |
| 13 | Compute normalization отсутствует | Центральный аргумент DART-Eval: «no compelling gains **per compute**». Без сравнения compute (GENA-LM pre-training ≈ $100–300k vs HyenaDNA ≈ $1–5k) claim про «overparameterization» не нормализован | **Medium** | Добавить одну строку в discussion: «GENA-LM pre-training compute ≈ 100× HyenaDNA; the 277× parameter ratio does not capture the full asymmetry» |
| 14 | E1 best layer per task — claim «иерархия features» не подтверждается | Если иерархия (k-mer → motif → regulatory), то peak должен сдвигаться вверх с complexity. Реально: peak Coding=5, Promoter=4, Enhancer=1, OCR=2 — это не монотонная иерархия | **Medium** | Не утверждать «hierarchy of features». Достаточно: «task-specific layer dependence varies» |
| 15 | Statistical significance не оценена нигде | E4 Silhouette 0.164 vs 0.124 = разница 0.04. Bootstrap CI? Permutation test? Без этого все «differences» необоснованы | **Medium** | Добавить bootstrap CI (200 resamples, 30 sec each) к E4 silhouette и E3 F1 |
| 16 | Отсутствуют ключевые citations (HyenaDNA, EPDnew, JASPAR, DeepSEA, BPE) | Academic correctness; ревьюер может счесть это небрежностью | **Medium** | Добавить 5 ссылок (см. раздел 4.6) |
| 17 | Multiple comparisons без correction | 4 tasks × 13 layers × 3 seeds × несколько метрик = много тестов. P(at least one «significant» finding by chance) высока | **Low** | Хотя бы упомянуть в limitations |

---

## 4. Detailed Section-by-Section Review

### 4.1 Title

**Текущий:** *GENA-LM: A Critical Review with Layer-Aware Reproducibility Experiments*

**Что работает:**
- Сразу понятно что это критика, не оригинальная статья
- «Layer-aware» — точный технический термин
- «Reproducibility» сигнализирует серьёзность

**Что не работает:**
- «Reproducibility» вводит в заблуждение — proposal **не** воспроизводит результаты статьи. Воспроизведение = repeat the same numbers. Proposal делает **extensions** (new probing analyses, new comparisons). Это **independent evaluation / extension**, не reproducibility.
- Title не передаёт самого сильного результата — HyenaDNA-vs-GENA-LM head-to-head.

**Альтернативные варианты:**
1. *GENA-LM: A Layer-Aware Critique with Frozen-Probing and Head-to-Head Evaluation* — точнее
2. *Where Does GENA-LM Encode Biology? A Layer-Wise Probing Study with HyenaDNA Comparison* — research-question-driven
3. *Frozen-Probing Critique of GENA-LM: Layer Specialization, Counterfactual Motif Probing, and DART-Eval-Inspired Baseline Calibration* — полное описание содержимого
4. *Layer-Aware Critique and Lightweight-Baseline Calibration of GENA-LM* — компактно

**Рекомендация:** вариант (2) — самый читабельный и сразу задаёт исследовательский вопрос.

### 4.2 Introduction / Background (Section 1)

**Сильные стороны:**
- Чёткая литературная контекстуализация: DNABERT (512) → DNABERT-2 (4k) → NT v2 (12k) → HyenaDNA (1M) → GENA-LM. Линейный нарратив.
- Сразу указано venue и DOI — академическая аккуратность.
- Указание на CORE A* контекст (Caduceus, DART-Eval) даёт ревьюеру понять что автор знает поле.

**Слабые стороны:**
- Утверждение «HyenaDNA scales to 1 Mb but loses accuracy on benchmarks» — без citation. Это claim, требующий ссылки на конкретные числа.
- Нет явной motivation **зачем эта школа** должна заинтересоваться этой работой. Связь с лекциями школы (Кардымон, Бурцев — соавторы GENA-LM) **не упомянута** в Introduction, хотя это прямой связующий аргумент к мотивации поступления.

**Missing context:**
- Один абзац про то, почему именно AIRI bioinformatics группа должна оценивать эту работу: Кардымон руководитель группы, GENA-LM — её flagship; кандидат демонстрирует прямую релевантность профилю группы.

**Как переписать:**
Добавить 2–3 строки после venue note:
> «I selected GENA-LM specifically because (i) it is the flagship DNA foundation model from AIRI's own Bioinformatics group, making this proposal directly aligned with the research focus of the host institution; (ii) it occupies a unique position as the only open-source family combining BPE, sparse attention, and RMT; (iii) its design choices (BPE over single-nucleotide, transformer over SSM) sit on the central axis of current debate in DNA LM literature.»

### 4.3 Research Problem

**Замечание о формате:** в paper-critique proposal классическое «Research Problem» заменяется на «What's missing/wrong in the reviewed paper». Это **раздел 4 (Weaknesses)** в текущем proposal.

**Насколько проблема ясна:** Weaknesses table формулирует 5 конкретных проблем (W3, W6, W7, W10, W12) с numeric labels из ревью статьи. Хорошо.

**Где логические дыры:**
- W6 (cross-species drop) включён в Weaknesses table, но соответствующего эксперимента **нет** в Results. Раздел Weaknesses обещает 5 проблем, а решает 4 (W3 не закрывается напрямую, W6 в follow-ups). Несогласованность.
- W12 (no DART-Eval) обещано закрыть E6, но E6 не есть DART-Eval (см. Major Issue #1).

**Как усилить problem statement:**

Переписать вступление к разделу 4 (Weaknesses):
> «We identify five concrete weaknesses in GENA-LM as published; **the present proposal directly addresses W7, W10, and (partially) W12** through new experiments; W3 and W6 are discussed as natural follow-ups for the school project itself.»

Это **честно** ограничивает scope и показывает что автор знает разницу между «hand-waving коверкой» и «реально решённой проблемой».

### 4.4 Research Gap

**Есть ли настоящий gap:**
- W7 («no comparison with Caduceus/Evo/HyenaDNA») — **настоящий gap**, легко проверяется.
- W10 («token importance ≠ causality») — **настоящий gap**, авторы статьи сами это acknowledge.
- W12 («no DART-Eval evaluation») — gap **формально настоящий** (статья вышла раньше DART-Eval), но менее «вина авторов» — это temporal gap.
- Дополнительный gap (BERTology layer-wise probing) — **новый**, не упомянут в W-таблице, появляется только в разделе 5. Это правильный gap, но он **должен быть в W-таблице** (как W13, например).

**Доказан ли литературой:**
- Каждый W привязан к citation (DART-Eval, Caduceus, etc.). Хорошо.
- Но gap «no layer-wise BERTology» **не привязан** к existing BERTology literature (Tenney 2019, Vig 2021 — есть в References, но без явной связи с gap).

**Как сформулировать точнее:**
Добавить в Weaknesses table:

| W13 | No layer-wise analysis: authors report only species-level layer structure (Fig 4), not task-wise/functional layer-wise probing | Standard BERTology methodology (Tenney 2019, Vig 2021) is systematically missing for DNA LMs |

### 4.5 Research Questions / Objectives / Hypotheses

**Замечание:** в proposal **нет явных research questions**. Эксперименты E1–E6 описаны как «проверим X», но без явной формулировки RQ.

**Это слабость.** Reviewers ожидают увидеть 2–4 RQ перед методологией.

**Предложенные RQ (добавить в новый раздел 5.0 перед списком экспериментов):**

> **Research Questions**
> 1. **RQ1 (Layer specialization).** Does GENA-LM exhibit BERTology-style layer specialization for biological tasks, and if so, at which depth do features peak?
> 2. **RQ2 (Architectural alternatives).** In the frozen-probing regime on short-range regulatory tasks, does a 270×-smaller single-nucleotide model (HyenaDNA-tiny) achieve comparable or better performance than GENA-LM-base?
> 3. **RQ3 (Counterfactual motif grounding).** Has GENA-LM internalized the canonical CTCF binding grammar in a counterfactual (not merely correlational) sense?
> 4. **RQ4 (Baseline calibration).** On a regulatory-DNA classification proxy, do GENA-LM's frozen representations offer compelling gains over lightweight supervised baselines (k-mer + LogReg, TinyCNN)?

**Соответствие experiments:** RQ1 ↔ E1+E4, RQ2 ↔ E3, RQ3 ↔ E5, RQ4 ↔ E6. Полное mapping.

**Hypotheses (опционально, для строгости):**
- **H1** (под RQ1): peak F1 будет в mid-layer (4–8) для regulatory tasks, last layer ниже из-за MLM head alignment.
- **H2** (под RQ2): HyenaDNA-tiny ≥ GENA-LM-base на short-range tasks, где BPE низкое разрешение мешает.
- **H3** (под RQ3): cosine distance в embedding space будет коррелировать с PWM-predicted off-consensus penalty (r > 0.5), а **mutations вне мотива не дадут такой корреляции** (negative control).
- **H4** (под RQ4): GENA-LM frozen ≤ TinyCNN supervised on regulatory tasks по AUROC.

Pre-stated hypotheses **автоматически защищают от post-hoc cherry-picking critique**.

### 4.6 Literature Review

**Замечание о формате:** в 2-страничном AIRI proposal full literature review не ожидается; нужны 6–10 ключевых ссылок. В текущей версии 8 references — корректно по объёму.

**Что есть:**
- GENA-LM (target paper) ✅
- Caduceus, DART-Eval (CORE A* контекст) ✅
- Tenney 2019, Vig 2021 (BERTology methodology) ✅
- RMT (Bulatov 2022) ✅
- Evo 2 ✅
- Genomic Benchmarks ✅

**Что критически отсутствует:**
1. **HyenaDNA paper** (Nguyen et al., NeurIPS 2023) — **главный baseline в E3 и E6**, отсутствие citation — серьёзная academic ошибка.
2. **DeepSEA** (Zhou & Troyanskaya, *Nat Methods* 2015) — упоминается в Method, нет в References.
3. **JASPAR / HOCOMOCO** — в E5 используется MA0139.1 CTCF motif, источник PWM должен быть процитирован (JASPAR Castro-Mondragon et al., *NAR* 2022/2024 update).
4. **EPDnew** (Dreos et al., *NAR* 2017) — источник promoter data, должен быть в References.
5. **BPE** (Sennrich et al., ACL 2016) — упомянуто в method, нет citation.

**Как исправить:** Добавить 5 references выше. Объём 2 страниц допускает компактный bibliography style.

### 4.7 Theoretical / Conceptual Framework

**В proposal не выделено как отдельный раздел.** Имплицитный framework — BERTology + DART-Eval critique.

**Что нужно уточнить:**
Один абзац в начале раздела 5 (Proposed improvement):
> «Conceptual framework. I approach GENA-LM through two complementary lenses: (i) **BERTology** (Tenney et al. 2019; Vig et al. 2021) — the methodology of layer-wise probing developed for NLP transformers, asking *where* in a frozen model task-specific features accumulate; (ii) **DART-Eval critique** (Patel et al. 2024) — the empirical claim that DNA LMs require fair comparison with lightweight supervised baselines and architectural alternatives to justify their compute. Both frameworks are systematically missing from the original GENA-LM evaluation; my experiments apply both to the frozen `gena-lm-bert-base-t2t` checkpoint.»

Это превращает «список 5 экспериментов» в **theoretically grounded research programme**, что сильно усиливает proposal.

### 4.8 Methodology

**Соответствует ли метод RQ:** Frozen probing — стандартный method для BERTology RQ1, и это правильный выбор. Но:

**Methodological gaps:**

1. **Frozen probing only** — RQ2 (HyenaDNA vs GENA-LM) формально требует fine-tuned comparison для apple-to-apples сравнения с published numbers. Этого нет.
2. **No control для E5** — нет run с random-init GENA-LM (которая должна дать r ≈ 0). Без этого нельзя сказать что high r — это «выученная грамматика», а не artefact инициализации/токенизации.
3. **Pooling strategy не varied** — mean-pool используется для всех экспериментов, ablation нет.
4. **Sequence length не варьируется** — все эксперименты на default context, нет проверки sensitivity к 1k vs 4k vs 32k для GENA-LM.
5. **No hyperparameter search для LogReg** — `C=1, max_iter=2000` — fixed. На разных задачах optimal C может быть разным.

**Methodological strengths:**
- StandardScaler + LogReg — корректный, простой, reproducible probing.
- 70/30 stratified split — стандартный.
- 3 seeds для E1 — корректно.

**Что не описано в Experimental Setup table:**
- HyenaDNA version (HF ID)
- Точная sequence length для каждой модели в E3 и E6
- Какой слой использован для HyenaDNA (только last? best per task?)
- Pooling для HyenaDNA — mean? CLS? Last hidden state?

**Как усилить validity/reliability:**

| Усиление | Усилие | Эффект |
|---|---|---|
| 3 seeds для E3, E6 | 2 ч | Закрывает stat-significance critique |
| Bootstrap CI для E4 silhouette | 30 мин | Закрывает «is 0.164 vs 0.124 meaningful?» |
| Random-init GENA-LM control в E5 | 30 мин | Закрывает «is r=0.893 artefact?» |
| Pooling ablation (CLS/mean/max) на одной задаче | 1 ч | Закрывает «BPE vs single-nuc pooling asymmetry» |
| Fine-tuned GENA-LM на промоутере | 30 мин | Закрывает «frozen vs fine-tuned» critique |

**Ethical considerations:** для compute-only experiments с public data из HuggingFace ethics блок не нужен. Можно добавить одну строку: «All data used (Genomic Benchmarks, HuggingFace models) is publicly available under permissive licences; no human subjects involved.»

### 4.9 Data / Sample / Materials

**Достаточно ли описаны:**
- 4 датасета из Genomic Benchmarks указаны явно — ✅
- Размеры выборок: **отсутствуют**. Сколько samples per task? Стандартные splits Genomic Benchmarks? — нужно явно указать.
- E5 «100 synthetic sequences with CTCF motif» — указано в Results, но не в Setup table ⚠️
- E6 «5,000 sequences, 60/20/20 split» — указано в Results, но не в Setup table ⚠️

**Sample size discussion отсутствует.** Достаточен ли 5000 для AUROC differences порядка 0.03? Standard error of AUROC при n=1000 (test set) ≈ 0.01–0.015. То есть разница 0.030 = ~2 σ — статистически граничит со significant.

**Как уточнить:**

Добавить в Experimental setup (заменить нынешнюю таблицу или дополнить):

| Task | Train | Val | Test | Source |
|---|---:|---:|---:|---|
| Promoter (nontata) | XXX | XXX | XXX | Genomic Benchmarks `human_nontata_promoters` |
| Enhancer (Cohn) | XXX | XXX | XXX | Genomic Benchmarks `human_enhancers_cohn` |
| OCR (Ensembl) | XXX | XXX | XXX | Genomic Benchmarks `human_ocr_ensembl` |
| Coding vs intergenic | XXX | XXX | XXX | Genomic Benchmarks `demo_coding_vs_intergenomic_seqs` |
| E5 CTCF | 100 synthetic | — | — | JASPAR MA0139.1 + random flanks |
| E6 OCR (DART proxy) | 3000 | 1000 | 1000 | `human_ocr_ensembl` |

XXX — заполнить из скриптов.

### 4.10 Data Analysis Plan

**Что есть:**
- E1: LogReg per layer per task, метрики F1/MCC/AUC ✅
- E4: KMeans + Silhouette/NMI/ARI ✅
- E3: same as E1, two models ✅
- E5: cosine distance correlation with PWM ✅
- E6: 4 models, AUROC/AUPRC/F1/MCC ✅

**Что нужно добавить:**
1. **Statistical test plan** для главных сравнений:
   - E3: paired t-test или Wilcoxon signed-rank на 3 seeds × 4 tasks между моделями
   - E5: bootstrap CI для r (200 resamples)
   - E6: DeLong test для AUROC differences
2. **Multiple comparison correction** plan (Bonferroni или Benjamini-Hochberg для 4–8 main comparisons)
3. **Sensitivity analyses**:
   - E1: F1 stability к выбору LogReg C (try C ∈ {0.1, 1, 10})
   - E3: stability к pooling (CLS, mean, max)

### 4.11 Expected Contribution / Significance

**Текущая формулировка (Section 8):**
6 пунктов (i)–(vi). Это много для 2-страничного proposal и выглядит как inflation.

**Какие сильны:**
- (i) reproducible BERTology pipeline — конкретно
- (ii) first systematic layer-aware decomposition — overclaim («first»), нужно «to our knowledge»
- (iv) counterfactual evidence CTCF — strong, но требует negative control
- (v) 65-feature beats 110M model — strong sound-bite

**Какие слабы:**
- (iii) «empirical evidence that a 254×–1000× smaller model matches or beats GENA-LM» — overclaim в frozen-only setting
- (vi) practical guidelines — generic

**Не завышен ли вклад:**
- «First systematic layer-aware decomposition of GENA-LM» — нужно проверить (есть ли blog posts? GitHub issues?), либо смягчить до «To our knowledge, the first published systematic layer-aware analysis…»
- «Causal evidence that GENA-LM has internalized CTCF» — overclaim без контролей.

**Theoretical vs empirical vs methodological vs practical:** в текущем виде не разделены. Лучше структурировать:

> **Contributions.**
> *Methodological:* a reusable BERTology + frozen-probing pipeline for DNA foundation models, open-sourced (E1, E4 code).
> *Empirical:* (i) layer-wise specialization of GENA-LM features; (ii) head-to-head frozen comparison with HyenaDNA showing competitive performance of a 270×-smaller architecture on short-range regulatory tasks; (iii) DART-Eval-inspired baseline calibration replicating the «no compelling gains» finding on a regulatory-DNA proxy.
> *Interpretive:* (iv) saturation-mutagenesis evidence consistent with internalization of the CTCF binding consensus by GENA-LM, supported by a planned random-position control in follow-up.
> *Practical:* (v) recommendation to prefer mid-layer (not last layer) for unsupervised downstream use of GENA-LM, with the caveat that this preference holds for geometric arrangement metrics specifically.

### 4.12 Limitations

**В текущем proposal limitations explicitly не выделены.** Это серьёзный пропуск. Любой proposal должен иметь Limitations block.

**Минимальный limitations block (добавить как раздел 8 или 8.5):**

> **Limitations.**
> (i) All experiments use **frozen embeddings**; fine-tuned performance may differ substantially, and findings should not be over-extrapolated to fine-tuned settings where authors report stronger GENA-LM performance (e.g., F1=94.16 on promoter 2 kb).
> (ii) E6 is a **DART-Eval-inspired proxy benchmark**, not the formal DART-Eval Task 1–5 pipeline; full DART-Eval evaluation is in follow-up.
> (iii) Experiments are limited to `gena-lm-bert-base-t2t` (110 M, 4.5 kb context); behavior of the BigBird (36 kb context) and RMT (Mb context) variants — where GENA-LM has architectural advantage — is not tested.
> (iv) HyenaDNA comparison uses the smallest `hyenadna-tiny-1k-seqlen` variant; larger HyenaDNA variants may compress the gap differently.
> (v) Sample sizes (≈2,000–5,000 sequences per task) limit statistical power; reported metric differences in E3 (Δ F1 ≈ 0.01–0.05) are at the edge of seed-noise resolution.
> (vi) All four tasks in E3 are short-range regulatory; long-range tasks (species classification 32 kb, DeepSEA HM 8 kb), where GENA-LM has documented advantage, are not tested in this proposal.

**Это критически важно:** ревьюер AIRI **точно** ожидает увидеть Limitations. Их отсутствие читается как «автор не понимает границ своей работы».

### 4.13 Timeline / Feasibility

**В proposal timeline не явен.** Сказано «≈ 16 minutes total for all five experiments» в Contribution — это compute time, не project timeline для самой школы.

**Что нужно добавить (если требуется AIRI форматом):**

> **Timeline for summer school (if accepted).**
> Week 1: replicate baseline GENA-LM probing pipeline on cluster; add `bert-large-t2t` (24 layers, 336 M); add `bigbird-base-t2t` (36 kb).
> Week 2: full DART-Eval Task 1–2 evaluation; Caduceus baseline (once mamba_ssm compiles on V100).
> Week 3: cross-species few-shot LoRA experiments (closing W6); E5 additional motifs (GATA2, ATF1).
> Week 4: writeup, ablations, supplementary experiments.

**Feasibility:** компактный proposal с realistic compute (V100 × 16 минут) — feasibility сильна. Это плюс.

---

## 5. Hidden Gaps and Blind Spots

Эти проблемы могут быть **неочевидны автору**, но опытный ревьюер их найдёт.

### 5.1 Скрытая assumption: «frozen ≈ fine-tuned для сравнительных claims»

Весь central narrative «GENA-LM не лучше HyenaDNA» неявно предполагает, что **frozen relative ordering моделей предсказывает fine-tuned ordering**. Это **не доказано** и часто **неверно**: некоторые модели «раскрываются» только при fine-tuning (особенно BPE-токенизаторы, где gradient flow через embeddings важен).

- **Где проявляется:** Section 7 E3 и E6 интерпретации.
- **Почему опасно:** ревьюер AIRI Bioinformatics — соавторы GENA-LM — точно знают что fine-tuned GENA-LM лучше frozen.
- **Как исправить:** явный disclaimer + одна fine-tuned точка для калибровки.

### 5.2 Conceptual confusion: «causal» vs «counterfactual» vs «correlational»

В E5 используется слово «causal». В строгом смысле causal claim требует **interventional manipulation на real системе**. То, что делается в E5 — это **counterfactual mutation на synthetic sequences**, что **сильнее correlational** (Figure 2 IG в оригинале), но **слабее true causal** (которое было бы fine-tuned model + real binding sites + functional perturbation измеренная в репортёрной системе).

- **Где проявляется:** Results E5, Contribution (iv).
- **Почему опасно:** академическая precision важна; «causal» — overclaim, который методолог поймает за полминуты.
- **Как исправить:** заменить «causal» на «counterfactual» или «interventional in silico» везде в E5.

### 5.3 Scope creep: 5 экспериментов в 2-страничном proposal

5 экспериментов (E1, E3, E4, E5, E6) — это много для 2 страниц. Каждый получает 100–150 слов в Results, что **не позволяет глубокого discussion**. Это создаёт ощущение «breadth-first, depth shallow».

- **Альтернатива:** **3 эксперимента глубоко** (E1, E3, E6 или E1, E3, E5) с детальной аргументацией, ablations, и controls для каждого.
- **Trade-off:** 5 экспериментов → больше «брэйкеров точек» в narrative (выглядит ambitious), но каждая слабее. 3 эксперимента → каждый solid, но менее ambitious.
- **Рекомендация:** для AIRI школы 3 глубоких лучше 5 поверхностных — школа отбирает по methodological rigor, не по count.

### 5.4 Confirmation bias в выборе задач для E3

Все 4 задачи в E3 — **short-range regulatory**. Это именно та территория, где BPE GENA-LM **должна** проиграть single-nuc HyenaDNA по architectural reasons (W3). Не включён ни один **long-range** task, где GENA-LM должна выиграть (species classification 32 kb, DeepSEA HM 8 kb).

- **Где проявляется:** E3 task selection.
- **Почему опасно:** ревьюер скажет «вы тестировали GENA-LM на её слабой территории и удивлены что она проиграла».
- **Как исправить:** добавить хотя бы одну long-range задачу. Либо явное disclaimer в Limitations: «E3 tasks are deliberately short-range to isolate the single-nucleotide-resolution effect; long-range performance, where GENA-LM has architectural advantage, is not tested in the proposal.» (Это уже включено в предложенный Limitations block.)

### 5.5 Measurement problem: cosine distance в E5 и BPE токенизация

Cosine distance в high-dim embedding space (768-dim для bert-base) **может быть несимметрична**: малые изменения нуклеотида при BPE могут приводить к **discontinuous** изменениям токенизации (одна замена → меняется граница BPE токенов → меняется hidden representation целого фрагмента). Это **не «модель не понимает CTCF»**, это **artefact BPE токенизации**.

- **Где проявляется:** E5 r=0.893.
- **Почему опасно:** the high r может быть **полностью** артефактом BPE, не learned biology.
- **Как исправить:**
  1. Negative control: same procedure с random-init модели → должно дать r ≈ 0, если r=0.893 — learned signal.
  2. Доп. check: для каждого мутанта проверить **изменилась ли токенизация**. Если каждая off-consensus мутация даёт BPE-shift, а каждая consensus не даёт — r=0.893 артефакт.
  3. Доп. control: мутации в **случайной позиции вне мотива** → должны дать r ≈ 0.

### 5.6 Слабая связь между literature и methodology

Tenney 2019 и Vig 2021 (BERTology references) **названы**, но их **methodology детально не применена**. Tenney использует **edge probing**, не linear probing. Vig анализирует **attention patterns**, не hidden states. Линковка к classic BERTology methodology формальная, не essential.

- **Где проявляется:** Section 5 (Proposed improvement), References.
- **Почему опасно:** ревьюер скажет «вы цитируете BERTology, но не используете её методы».
- **Как исправить:** либо добавить attention-pattern analysis (соответствие Vig), либо признать что это «BERTology-inspired», not «BERTology methodology proper». Минимум — сказать честно: «We adopt linear-probing as a tractable subset of the broader BERTology toolkit; attention-pattern analysis (Vig 2021) is left for follow-up.»

### 5.7 Риск: что вообще «meaningful» finding значит здесь?

Все 5 экспериментов дают численные результаты, но **что было бы negative finding**? Если бы layer-wise probing показал плоскую кривую (нет peak) — это было бы negative finding или другое positive? Если бы HyenaDNA проиграла GENA-LM на всех 4 — это закрыло бы W7? Pre-registration of hypotheses отсутствует, что создаёт риск **post-hoc interpretation** — «whatever we found, we'll explain it».

- **Где проявляется:** Section 5 + Section 7.
- **Почему опасно:** все «findings» можно интерпретировать как «interesting», что подозрительно для ревьюера.
- **Как исправить:** добавить pre-stated hypotheses (см. раздел 4.5 этого ревью) — это автоматически делает findings interpretable как «hypothesis confirmed / rejected».

### 5.8 Operational definition of «mid-layer»

В тексте используется «mid-layer» как understood termin, но **operationally** mid-layer для bert-base (12 layers) — это L4-L8? L5-L7? L6? Точная operational definition отсутствует. В разных результатах «mid» означает разное (L4, L5, или просто «not L0 and not L12»).

- **Как исправить:** одна строка в Methodology: «Throughout, *mid-layer* refers to layers L4–L7 (one-third to two-thirds depth); *embedding layer* refers to L0; *last layer* refers to L12.»

---

## 6. Recommended Fixes

### 6.1 Immediate Fixes (4–5 часов работы, MAX impact)

| # | Fix | Время | Эффект |
|---|---|---|---|
| 1 | Перезапустить E3 и E6 с n=3 seeds, добавить std в таблицы | 2 ч | Закрывает Major Issue #2 — главный stat-significance критика |
| 2 | Переименовать E6 «DART-Eval» → «DART-Eval-inspired proxy» везде | 15 мин | Закрывает Major Issue #1 — главный overstating |
| 3 | Добавить «frozen-probing regime» qualifier в каждый claim про HyenaDNA winning (минимум 4 места: Abstract, Results E3/E6, Contribution) | 30 мин | Закрывает Major Issue #3 |
| 4 | E5 negative control: 200 forward passes с мутациями в случайных позициях flank вне мотива, добавить r-row в таблицу | 30 мин | Закрывает Major Issue #4 — causal/counterfactual claim под контролем |
| 5 | Уточнить HyenaDNA HF ID и точную sequence length для обеих моделей в Setup table | 15 мин | Закрывает Major Issue #10 |
| 6 | Добавить Limitations block (см. 4.12 этого ревью) | 30 мин | Закрывает Major Issue #8 — критический пропуск раздела |
| 7 | Заменить «causal» на «counterfactual» в E5 везде | 15 мин | Закрывает 5.2 overclaim |
| 8 | Добавить HyenaDNA, EPDnew, JASPAR, DeepSEA, BPE citations | 20 мин | Закрывает Major Issue #16 — academic correctness |

**Итого: ~4 часа 35 минут.**

### 6.2 Structural Fixes (если есть полдня — 4–6 часов дополнительно)

| # | Fix | Время | Эффект |
|---|---|---|---|
| 9 | Добавить явные Research Questions (4 RQ + 4 H) в начало Section 5 | 1 ч | Делает proposal hypothesis-driven, не результат-driven |
| 10 | Переписать E1 интерпретацию: «peak в mid-layer для задач что модель решает; near-chance для long-range» | 30 мин | Закрывает Major Issue #5 |
| 11 | E4 честное обсуждение всех трёх метрик (Silhouette/NMI/ARI), не cherry-pick | 20 мин | Закрывает Major Issue #6 |
| 12 | E3 дополнительно с GENA-LM на best-per-task layer (из E1) — таблица расширяется | 1 ч | Закрывает Major Issue #7 |
| 13 | Restructure Contributions section по methodological/empirical/interpretive/practical | 30 мин | Делает Contribution прозрачнее |
| 14 | Добавить Conceptual Framework абзац | 20 мин | Делает proposal theoretically grounded |
| 15 | Заполнить Sample sizes в Setup table | 10 мин | Закрывает «numbers from nowhere» риск |

**Итого структурного: +3.5 часа сверх immediate.**

### 6.3 Advanced Improvements (если есть полный день — ~8 часов сверх)

| # | Fix | Время | Эффект |
|---|---|---|---|
| 16 | Одна fine-tuned точка GENA-LM на промоутере (30 min compute) + сравнение с frozen | 1.5 ч | Колоссально усиливает proposal, закрывает frozen-vs-fine-tuned критику |
| 17 | E5 на random-init GENA-LM как контроль | 1 ч | Strengthens counterfactual claim до почти-causal |
| 18 | Pooling ablation: CLS vs mean vs max для GENA-LM и HyenaDNA на одной задаче | 1 ч | Закрывает Major Issue #11 |
| 19 | Statistical test plan: paired Wilcoxon для E3, DeLong для E6, bootstrap CI для E4 | 1 ч | Усиливает статистическую строгость |
| 20 | Добавить одну long-range задачу в E3 (например species classification 8 kb) | 2 ч | Закрывает confirmation bias критика (5.4) |
| 21 | bert-large probing (вместо bert-base) для E1 — bonus | 1.5 ч | Демонстрирует scaling — другой peak pattern? |

---

## 7. Revised Research Design

### 7.1 Refined Topic
*Layer-Aware Frozen-Probing Critique of GENA-LM with Architectural and Baseline Calibration*

### 7.2 Revised Problem Statement

> «GENA-LM (Fishman et al., NAR 2025) represents the current state-of-the-art open-source DNA foundation model family from AIRI, combining BPE tokenization, sparse attention, and Recurrent Memory Transformer to scale context to megabases. However, its evaluation in the original paper has three systematic gaps relevant to the BERTology + DART-Eval era: **(i)** no task-wise layer-wise analysis of where in the network biological information accumulates — a standard BERTology question that determines which layer downstream users should query; **(ii)** no head-to-head comparison with non-BPE single-nucleotide architectures (HyenaDNA, Caduceus, Evo) that have emerged as direct alternatives; **(iii)** no calibration against the DART-Eval (NeurIPS 2024) critique that current DNA LMs may not justify their compute over lightweight supervised baselines. The proposal addresses gaps (i)–(iii) through five reproducible frozen-probing experiments on a single Tesla V100, with explicit acknowledgement of the frozen-only scope.»

### 7.3 Revised Research Gap

> «To my knowledge, no published analysis of GENA-LM combines: (a) task-wise layer-wise BERTology methodology applied to DNA foundation models; (b) head-to-head frozen comparison with single-nucleotide alternatives on regulatory tasks; (c) DART-Eval-inspired baseline calibration with explicit comparison to lightweight supervised models. The proposed experiments fill this gap on a constrained (≤ 1 V100, ≤ 20 min compute) reproducible budget.»

### 7.4 Revised Research Questions

См. раздел 4.5 этого ревью (RQ1–RQ4 + H1–H4).

### 7.5 Revised Objectives

1. **O1:** Characterize layer-wise feature specialization in frozen GENA-LM across four binary regulatory tasks (E1).
2. **O2:** Test whether layer-wise patterns generalize to unsupervised geometry (E4).
3. **O3:** Benchmark GENA-LM-base frozen representations against HyenaDNA-tiny on identical tasks and protocol (E3).
4. **O4:** Test whether GENA-LM has internalized the CTCF binding consensus through counterfactual saturation mutagenesis with negative controls (E5).
5. **O5:** Calibrate GENA-LM against lightweight supervised baselines on a DART-Eval-inspired regulatory-DNA proxy task (E6).

### 7.6 Recommended Methodology

См. Section 4.8 выше — frozen probing + 3 seeds + statistical tests + controls. Plus один fine-tuned datapoint для калибровки (если есть полдня).

### 7.7 Expected Contribution

См. revised version в разделе 4.11 — 4 категории (methodological / empirical / interpretive / practical).

---

## 8. Suggested Rewrites

### 8.1 Improved Abstract / Opening

**Текущее (имплицитное в Section 1):** «Earlier DNA LMs are short-context, HyenaDNA loses accuracy at scale, GENA-LM closes the gap.»

**Предложенное (новый Abstract):**
> «GENA-LM (Fishman et al., NAR 2025) is the flagship open DNA foundation model from AIRI Bioinformatics, but its evaluation predates two methodological developments that fundamentally reshape how DNA LMs should be assessed: **layer-wise BERTology probing** and the **DART-Eval critique** of unjustified scale. In this proposal, I apply both frameworks to the frozen `gena-lm-bert-base-t2t` checkpoint across four regulatory-DNA classification tasks (`promoter`, `enhancer`, `OCR`, `coding`), one interventional motif-probing experiment on CTCF, and one DART-Eval-inspired baseline calibration. All experiments run in under 20 minutes on a single Tesla V100. **Key findings (frozen-probing regime).** (i) Mid-layer specialization for tasks the model meaningfully resolves (promoter, coding); long-range regulatory tasks (enhancer, OCR) remain near-chance at any layer, suggesting that BigBird and RMT variants are required. (ii) A 270×-smaller single-nucleotide alternative (HyenaDNA-tiny) matches or exceeds GENA-LM on three of four short-range regulatory tasks. (iii) A 65-feature 3-mer + LogReg baseline approaches GENA-LM (Δ AUROC = +0.011) on Open Chromatin classification — an empirical replication of the DART-Eval critique on a regulatory-DNA proxy. (iv) Counterfactual saturation mutagenesis on CTCF reveals embedding responses consistent with the canonical binding consensus (r = 0.893 against PWM-derived expectations), pending negative-control validation. Limitations and follow-ups for the summer school are discussed.»

### 8.2 Improved Research Gap statement

> «Three concrete gaps motivate this proposal: **(G1)** the GENA-LM paper analyses layer structure only at the species-level (Fig 4) and does not perform task-wise layer-wise probing in the BERTology tradition (Tenney et al. 2019; Vig et al. 2021), leaving downstream users without guidance on which layer to use for unsupervised tasks; **(G2)** the paper compares GENA-LM with prior BPE/k-mer DNA LMs (DNABERT, NT v2, HyenaDNA on species classification only), but does not include a frozen-probing head-to-head against single-nucleotide alternatives (HyenaDNA on regulatory tasks, Caduceus, Evo 2) which form the central architectural axis of current debate; **(G3)** the paper predates the DART-Eval benchmark (Patel et al., NeurIPS 2024 D&B) and its central critique that DNA LMs may not justify their compute over lightweight supervised baselines.»

### 8.3 Improved methodology justification

> «I adopt a **frozen-probing** methodology for three reasons: **(i)** it isolates the representational content of pre-trained weights from task-specific gradient adaptation, which is the methodologically correct test for «what has the foundation model learned»; **(ii)** it matches the BERTology tradition (Tenney 2019, Vig 2021), which is the explicit methodological reference for RQ1; **(iii)** it is computationally tractable within a 1×V100 / 2-week budget. I explicitly acknowledge that fine-tuned performance may differ substantially (see Limitations) and provide one fine-tuned datapoint for calibration on the promoter task.»

### 8.4 Improved aim statement

> «**Aim.** To produce a layer-aware, baseline-calibrated, counterfactually-probed characterization of GENA-LM's frozen representational content, situated within the BERTology and DART-Eval literatures, and to use the resulting empirical picture to identify which architectural and methodological extensions are most promising for follow-up during the summer school project itself.»

### 8.5 Improved contribution statement

См. раздел 4.11 — restructured по methodological/empirical/interpretive/practical.

### 8.6 Improved limitations statement

См. раздел 4.12 — полный 6-point block.

### 8.7 Improved «causal CTCF» phrasing (E5)

**Текущее:**
> «**Causal** interpretation counterpart to Figure 2 of the paper — directly closes W10.»

**Предложенное:**
> «A **counterfactual / in silico interventional** complement to the **correlational** integrated-gradients analysis in Figure 2. While not causal in the strict experimental sense (which would require functional perturbation in a reporter assay), it represents a stronger interpretive standard than attribution methods: each substitution is an actual change to the input, and the embedding response is measured directly. Combined with the planned negative control (mutations at random positions outside the motif, expected r ≈ 0), this supports the claim that GENA-LM has internalized the CTCF binding consensus rather than reacting to BPE-tokenization artefacts.»

### 8.8 Improved «HyenaDNA beats GENA-LM» phrasing

**Текущее:**
> «HyenaDNA wins 3/4 tasks with 277× fewer parameters»

**Предложенное:**
> «**In the frozen-probing regime, on these four short-range regulatory tasks**, HyenaDNA-tiny achieves comparable-or-better F1 than GENA-LM-base on three out of four tasks (promoter, enhancer, OCR), while using approximately 277× fewer parameters and an order of magnitude less pre-training compute. GENA-LM-base leads on the coding-vs-intergenic task, consistent with the strength of BPE / k-mer compositional features there. This pattern — small specialized models matching larger foundation models on short-range regulatory tasks — directly replicates the spirit of the DART-Eval critique (Patel et al., NeurIPS 2024). It does not generalize to long-range tasks (e.g., species classification at 32 kb), where GENA-LM has documented architectural advantage and is not tested here.»

---

## 9. Questions I Must Answer Before Final Submission

### Conceptual clarity
- **Q1.** Что именно ты defends если ревьюер скажет «frozen vs fine-tuned — это не fair сравнение моделей»? Готовь explicit response (1–2 sentences).
- **Q2.** Что является **negative finding** для каждой из 5 экспериментов? Если не можешь сформулировать, hypotheses не pre-stated.
- **Q3.** Почему именно эти 4 задачи, а не другие из Genomic Benchmarks? (есть `human_ensembl_regulatory`, `drosophila_enhancers_stark`, `demo_human_or_worm` — почему не они?)
- **Q4.** Operational definition of «mid-layer» в твоей работе — это L4 specifically? L4–L7 range? Уточнить.

### Literature
- **Q5.** Есть ли HyenaDNA paper в References? (нет — добавить Nguyen et al. NeurIPS 2023)
- **Q6.** Есть ли EPDnew, JASPAR, DeepSEA, BPE citations? (вероятно нет — добавить)
- **Q7.** Цитировал ли ты *официальный* Genomic Benchmarks paper (Grešová) для каждой задачи? Особенно для нестандартного `human_nontata_promoters` subset?

### Methodology
- **Q8.** Сколько seeds в E3 и E6? Если 1 — перезапустить.
- **Q9.** Какой именно `hyenadna-tiny-*` checkpoint используется? Какая max sequence length для каждой модели?
- **Q10.** Mean-pooling делается **до или после** masking padding tokens? Это влияет на BPE asymmetry.
- **Q11.** В E5 — что происходит с BPE токенизацией при single-nucleotide substitution? Изменяется ли разбиение на токены? (можно проверить grep'ом тестового скрипта).
- **Q12.** Какой слой использован для HyenaDNA в E3/E6 — last? best-per-task? (для GENA-LM используется last, что является weak choice — см. Issue #7).

### Data
- **Q13.** Точные train/val/test размеры для каждой из 4 задач?
- **Q14.** В E5 — flanks случайные (random ACGT) или из real intergenic регионов? (это влияет на «causality» framing).
- **Q15.** В E6 — 60/20/20 split стратифицирован?

### Feasibility (для самой школы)
- **Q16.** Есть ли realistic plan для DART-Eval Task 1–5 если accepted? Нужен Synapse access (~1 неделя bureaucracy).
- **Q17.** Что делаешь если mamba_ssm всё ещё не компилируется на V100 к началу школы? Backup plan для Caduceus.

### Contribution
- **Q18.** Что в твоей работе **first** vs «to my knowledge first» vs derivative от existing approaches? Honest accounting перед submission.
- **Q19.** Если ревьюер скажет «это просто replication DART-Eval в miniature» — что отличает твою работу? (Подсказка: layer-aware angle + counterfactual CTCF + AIRI-specific model focus.)

### Ethics
- **Q20.** Использовал ли ты any human-derived data, кроме reference genome / chromatin tracks? Для всех Genomic Benchmarks tasks ответ — нет, всё public reference data, но явно сказать в Methodology.

---

## 10. Final Verdict

### Можно ли подавать в текущем виде

**Да, можно — но ниже потенциала.** Сейчас proposal **пройдёт первый фильтр**, но **не максимизирует балл** от компетентного ревьюера. Ревьюером с высокой вероятностью будет сотрудник AIRI Bioinformatics — то есть соавтор GENA-LM или ближайшая группа. Они **точно** заметят:
- frozen vs fine-tuned confusion;
- DART-Eval подмена (Genomic Benchmarks ≠ DART-Eval);
- отсутствие seeds для E3/E6;
- overclaim про causal CTCF без контролей;
- отсутствие Limitations block.

Каждое из этих — это **conversation starter на интервью** (если будет), либо **silent downgrade** в оценке (если бумажный отбор).

### Главные риски при оценке

1. **Риск «overstating»:** 3–4 центральных claims сформулированы строже, чем подтверждено данными. Ревьюер с domain expertise это увидит за 5 минут.
2. **Риск «не понимает scope»:** frozen probing интерпретируется как universal model comparison.
3. **Риск «cherry-picking»:** все 4 задачи в E3 — short-range, что выгодно для HyenaDNA. Не сбалансировано long-range задачей.
4. **Риск «no controls»:** E5 causal claim без negative control воспринимается как методологическая наивность.
5. **Риск «нет Limitations»:** прямой пропуск ожидаемого раздела.

### 5 исправлений с максимальным эффектом

1. **Перезапустить E3+E6 с 3 seeds** (2 ч) — закрывает stat-significance.
2. **«DART-Eval» → «DART-Eval-inspired»** + явный disclaimer (15 мин) — закрывает главный overstating.
3. **Frozen-probing qualifier везде** (30 мин) — закрывает apples-vs-oranges.
4. **E5 negative control** (30 мин) — закрывает causal overclaim.
5. **Limitations block** (30 мин) — закрывает критический пропуск.

**Итого: ~4 часа работы, и proposal переходит из «acceptable» в «strong».**

### Какой должна быть следующая версия

**v2 структура (рекомендация):**

1. Problem and selected paper (как сейчас, +абзац про AIRI alignment)
2. Method (как сейчас)
3. Strengths (как сейчас)
4. Weaknesses (расширить — добавить W13 «no layer-wise»)
5. **Research Questions + Hypotheses (новый раздел)** — 4 RQ + 4 H
6. **Conceptual Framework (новый абзац)** — BERTology + DART-Eval lenses
7. Proposed experiments (как сейчас, но «E6 DART-Eval-inspired»)
8. Experimental setup (расширить — HF IDs, sample sizes, pooling spec)
9. Results (как сейчас + std + bootstrap CI + E5 negative control row)
10. **Limitations (новый раздел)** — frozen-only, short-range only, sample size, scope
11. Expected contribution (restructure по 4 категориям)
12. Follow-ups (как сейчас, + DART-Eval formal как priority #1)
13. References (добавить HyenaDNA, EPDnew, JASPAR, DeepSEA, BPE)

Объём: остаться в 2 страницы (или 2.5 для AIRI free format) — все добавления компактны при аккуратном редактировании.

---

### Финальная фраза

Работа сильная для уровня bachelor-application, но **на 4–5 часах подготовки можно поднять её с «competent» до «conspicuously rigorous»**. Это разница между «accepted» и «accepted in the top quartile». Учитывая что ты уже сделал **самую тяжёлую часть** — реальные эксперименты с реальными числами — оставить proposal в текущем виде после такой работы будет academically wasteful.

**Recommended action:** выделить один день (~5 часов) на immediate fixes (раздел 6.1) + добавить Limitations и Research Questions (раздел 6.2 пункты 9, 14). Это превратит proposal в submission, которым можно гордиться.
