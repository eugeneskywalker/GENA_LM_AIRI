---
stylesheet:
  - ./invitation.css
pdf_options:
  format: A4
  margin: 0
  printBackground: true
---

<div class="eyebrow"><span class="accent">Лето с AIRI · 2026</span> &nbsp; / &nbsp; Research Proposal</div>

<h1 class="title">GENA-LM</h1>

<div class="subtitle">анализ и развитие</div>

<p class="deck">Открытое семейство ДНК-моделей от AIRI, опубликованное в Nucleic&nbsp;Acids&nbsp;Research в 2025 году.</p>

<div class="cta cta-top">
  <div class="cta-label">Полный разбор с интерактивными графиками</div>
  <a class="cta-url" href="https://eugeneskywalker.github.io/GENA_LM_AIRI/">eugeneskywalker.github.io/GENA_LM_AIRI/</a>
</div>

<hr class="rule" />

<div class="section-label">1 · Метод статьи</div>

<p class="prose">
GENA-LM — открытое семейство ДНК-моделей от AIRI, опубликованное в Nucleic&nbsp;Acids&nbsp;Research в 2025 году. До неё ДНК-модели видели только короткие участки: DNABERT — 512 нуклеотидов, Nucleotide Transformer — 12&nbsp;000. Регуляторные элементы, влияющие на гены с расстояния 100&nbsp;000 букв, в такое окно не попадают. GENA-LM работает с участками <strong>до 36&nbsp;000 нуклеотидов сразу</strong>, а с рекуррентной памятью — до миллионов.
</p>

<p class="prose">
Главная техническая идея — не одно изобретение, а пайплайн из четырёх готовых компонентов, который впервые работает для длинных ДНК-последовательностей. Новизна — в их комбинации.
</p>

<ul class="bullets compact">
  <li><strong>BPE-токенизация вместо k-mers.</strong> Словарь на 32&nbsp;000 «слов» переменной длины (медиана 9 нт). В стандартное окно влезает в 8–10 раз больше ДНК.</li>
  <li><strong>BERT-архитектура с masked language modeling.</strong> 15&nbsp;% токенов прячутся — модель угадывает их по соседям. Никакой биологической разметки, только сама ДНК.</li>
  <li><strong>Sparse attention из BigBird.</strong> Локальные + случайные + глобальные связи. Длина входа вырастает с 4.5 до 36 кб.</li>
  <li><strong>Recurrent Memory Transformer (RMT).</strong> Длинная последовательность режется на сегменты по 4.5 кб; 10 memory-токенов «впитывают» прочитанное и передаются дальше — до миллионов нуклеотидов.</li>
</ul>

<p class="prose">
Данные — сборка T2T-CHM13v2 (полный геном человека, включая центромеры и теломеры) с вариантами из 1000&nbsp;Genomes. Для мультивидовых моделей — геномы мыши, дрозофилы, нематоды, дрожжей и десятков других эукариот.
</p>

<div class="section-label">2 · Ключевые задачи</div>

<p class="prose">
Семейство из <strong>8 предобученных моделей</strong> — base (110 млн параметров) и large (336 млн), плюс варианты с BigBird и RMT. Все веса, код, Colab и веб-сервис <em>dnalm.airi.net</em> — в открытом доступе. Авторы оценивают модели на <strong>9 биологических задачах</strong>, например «это промотор?». В моих экспериментах ниже я использую 4 из них: промоторы, энхансеры, открытая ДНК, кодирующие участки.
</p>

<div class="two-col">
  <div class="col">
    <div class="section-label">3 · Сильные стороны</div>
    <ul class="bullets compact">
      <li><strong>Открытость по полной программе.</strong> 8 моделей на HuggingFace, код, Colab, веб-сервис.</li>
      <li><strong>Действительно длинный контекст.</strong> На момент выхода — самая длинноконтекстная ДНК-модель на трансформере.</li>
      <li><strong>Cross-species transfer.</strong> F1&nbsp;≈&nbsp;0.95 на близких млекопитающих без переобучения.</li>
      <li><strong>Учит биологию.</strong> Сама выделяет мотивы транскрипционных факторов — ATF1, GATA2, CTCF.</li>
    </ul>
  </div>
  <div class="col">
    <div class="section-label">4 · Слабые стороны</div>
    <ul class="bullets compact">
      <li><strong>Один токен ≈ 9 нуклеотидов.</strong> Точечные мутации почти не видны: AUC&nbsp;0.66 на ClinVar SNV.</li>
      <li><strong>Не везде побеждает.</strong> На промоторах в окне 300&nbsp;п.&nbsp;о. проигрывает DNABERT, на сплайсинге — SpliceAI.</li>
      <li><strong>Хрупкая инфраструктура.</strong> Sparse attention требует CUDA + Triton + DeepSpeed; ломается между версиями драйверов.</li>
      <li><strong>Нет сравнения с SSM.</strong> Caduceus и Evo появились почти одновременно — но их в работе нет.</li>
    </ul>
  </div>
</div>

<div class="page-break"></div>

<div class="section-label">5 · Минимальные эксперименты</div>

<p class="prose">
Авторы дообучали модель под каждую задачу и читали ответ только с последнего слоя. Я проверил, что находится внутри модели <em>до</em> дообучения и на других слоях.
</p>

<div class="exp-card">
  <div class="exp-tag">5.1 · что знает каждый слой</div>
  <div class="exp-title">Не читай последний слой.</div>
  <p class="exp-body">
  Снял эмбеддинги с каждого из 12 слоёв (плюс с входа L0), проверил на 4 задачах. <strong>Биология модели живёт в середине</strong> — около 4–5 слоя. Promoter: L4 = 0.809 vs L12 = 0.778. Coding: L5 = 0.915 vs L12 = 0.887.
  </p>
</div>

<div class="exp-card">
  <div class="exp-tag">5.2 · gena-lm vs hyenadna</div>
  <div class="exp-title">Маленькая модель «победила» большую. Спойлер — нет.</div>
  <p class="exp-body">
  GENA-LM против HyenaDNA, которая в ~280 раз меньше по числу параметров. С последнего слоя HyenaDNA выигрывает 3 из 4 задач. С лучшего слоя на каждую задачу — GENA-LM возвращает себе 3 победы. <strong>Вывод сравнения меняется в зависимости от выбора слоя.</strong>
  </p>
</div>

<div class="exp-card">
  <div class="exp-tag">5.3 · проверка причины</div>
  <div class="exp-title">Модель выучила биологию, а не запомнила «редкие буквы».</div>
  <p class="exp-body">
  Поменял буквы внутри мотива CTCF (биологически важный участок) — модель сильно реагирует, r&nbsp;=&nbsp;0.893 с PWM-важностью. В случайных местах вне мотива — почти не реагирует. <strong>Разница 15×.</strong> Это исключает гипотезу «реагирует на любые редкие токены».
  </p>
</div>

<hr class="rule rule-tight" />

<div class="section-label">6 · Что предлагаю развить</div>

<div class="roadmap">
  <div class="rm-pair">
    <div class="rm-find">
      <div class="rm-find-title">Лучший слой — не последний.</div>
      <p>GENA-LM из середины обходит HyenaDNA на 3 задачах из 4. С последнего слоя — проигрывает.</p>
    </div>
    <div class="rm-arrow">→</div>
    <div class="rm-prop">
      <div class="rm-prop-title">Сделать «лучший слой» стандартом.</div>
      <p>Измерять ДНК-модели по их лучшему слою, не последнему. Ожидаемый эффект: +3–5&nbsp;% точности.</p>
    </div>
  </div>

  <div class="rm-pair">
    <div class="rm-find">
      <div class="rm-find-title">Модель действительно знает биологию.</div>
      <p>15× разница между важными и случайными местами. Это не угадывание.</p>
    </div>
    <div class="rm-arrow">→</div>
    <div class="rm-prop">
      <div class="rm-prop-title">Повторить тест на других участках.</div>
      <p>Ещё 3 биологических сигнала и другие ДНК-модели. Общее ли это свойство — или повезло с CTCF.</p>
    </div>
  </div>
</div>

<hr class="rule rule-tight" />

<div class="section-label">Разбираемая статья</div>

<p class="citation">
Fishman V., Kuratov Y., Shmelev A., Petrov M., Penzar D., Shepelin D., Chekanov N., Kardymon O., Burtsev M.
<a class="paper-title" href="https://doi.org/10.1093/nar/gkae1310"><em>GENA-LM: a family of open-source foundational DNA language models for long sequences.</em></a>
Nucleic Acids Research <strong>53(2)</strong>, gkae1310 (2025).
<a href="https://doi.org/10.1093/nar/gkae1310" class="inline-doi">doi:10.1093/nar/gkae1310</a>
</p>

<div class="footer">
  Код: <a href="https://github.com/eugeneskywalker/GENA_LM_AIRI">github.com/eugeneskywalker/GENA_LM_AIRI</a>
</div>
