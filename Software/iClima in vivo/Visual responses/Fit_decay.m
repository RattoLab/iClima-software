%%Script fot the fit of the 5s  grid visual responses

% t=[];
% STEP=[];
% GRID5=[];
% GRID10=[];

%% fit_logistic_decay.m
% Fit tripartito con logistica in salita E in discesa
% per stimoli STEP, GRID5, GRID10.
%
% Variabili richieste nel workspace:
%   t       - vettore colonna [N x 1], tempo (s), uguale per tutte le matrici
%   STEP    - matrice [N x nAnimals]
%   GRID5   - matrice [N x nAnimals]
%   GRID10  - matrice [N x nAnimals]
%
% Modello tripartito:
%   t < 0            : M(t) = 0
%   0 <= t < t_peak  : M(t) = p1 * (1 - exp(-p2*t))^p3       [logistica salita]
%   t >= t_peak      : M(t) = A_peak / (1 + exp(k*(t-t_half)))[logistica discesa]
%
% Strategia fit (come nel codice precedente con decay esponenziale):
%   1. Smooth della traccia per trovare il picco in modo robusto
%   2. Fit separato salita  -> p1, p2, p3
%   3. Fit separato discesa -> k, t_half  (con A_peak fissato dal fit salita)
%   Questo e' molto piu' stabile di un fit simultaneo a 5 parametri.

%% ---- Configurazione ----
stimNames  = {'STEP', 'GRID5', 'GRID10'};
stimData   = {STEP,   GRID5,   GRID10};
nStim      = numel(stimNames);

% Finestra di smoothing per la ricerca del picco (in campioni)
smoothWin  = 0;   % regola se necessario

% Opzioni ottimizzatore
opts_lsq = optimoptions('lsqcurvefit', ...
    'Display',                'off', ...
    'MaxFunctionEvaluations', 5000, ...
    'MaxIterations',          2000, ...
    'FunctionTolerance',      1e-10, ...
    'StepTolerance',          1e-10);

% Colori animali
nAnimalsMax  = max(cellfun(@(x) size(x,2), stimData));
animalColors = lines(nAnimalsMax);

%% ---- Storage risultati ----
allR2     = cell(nStim, 1);
allParams = cell(nStim, 1);

%% ---- Loop stimoli ----
for iStim = 1:nStim

    data      = stimData{iStim};   % [nTime x nAnimals]
    nTime     = size(data, 1);
    nAnimals  = size(data, 2);
    stimLabel = stimNames{iStim};

    % Vettore tempo per questo stimolo
    % IMPORTANTE: usa i primi nTime campioni di t per evitare mismatch
    t_stim = t(1:nTime);

    % Indici baseline e post-stimolo
    baselineIdx = t_stim < 0;
    postIdx     = t_stim >= 0;
    t_post      = t_stim(postIdx);

    R2_stim     = nan(1, nAnimals);
    params_stim = cell(1, nAnimals);
    clear stimResults   % evita conflitti con iterazioni precedenti

    % ---------- Figura ----------
    nCols = ceil(sqrt(nAnimals));
    nRows = ceil(nAnimals / nCols);
    figH  = figure('Name', ['Fit - ' stimLabel], ...
                   'NumberTitle', 'off', ...
                   'Color', 'w', ...
                   'Units', 'normalized', ...
                   'OuterPosition', [0 0 1 1]);

    for iAnim = 1:nAnimals

        y_raw  = data(:, iAnim);

        % --- Sottrai baseline ---
        y_raw = y_raw - mean(y_raw(baselineIdx), 'omitnan');

        % --- Smooth SIMMETRICO solo per trovare il picco ---
        % movmean([k k]) usa k campioni prima e k dopo -> nessun ritardo di fase
        y_smooth = movmean(y_raw, [smoothWin smoothWin], 'omitnan');
        y_smooth = y_smooth - mean(y_smooth(baselineIdx), 'omitnan');

        % --- Trova il picco sulla traccia smoothata, solo post-stimolo ---
        y_post_smooth = y_smooth(postIdx);
        [peak_val_smooth, peak_idx_local] = max(y_post_smooth);
        t_peak = t_post(peak_idx_local);

        % Valore di picco dalla traccia raw (per il fit)
        y_post_raw = y_raw(postIdx);
        Apeak_raw  = y_post_raw(peak_idx_local);

        % --- Maschera salita e discesa ---
        riseMask  = postIdx & (t_stim <= t_peak);
        decayMask = t_stim >= t_peak;

        t_rise  = t_stim(riseMask);
        y_rise  = y_raw(riseMask);
        t_decay = t_stim(decayMask) - t_peak;   % scala da 0
        y_decay = y_raw(decayMask);

        % ================================================================
        % FIT 1: salita  ->  p1*(1-exp(-p2*t))^p3
        % ================================================================
        riseModel = @(p, t) p(1) .* (1 - exp(-p(2) .* t)) .^ p(3);

        p1_0 = Apeak_raw * 1.1;
        p2_0 = 3 / max(t_peak, 0.1);
        p3_0 = 2;
        p0_rise = [p1_0, p2_0, p3_0];

        % Bounds: p1 puo' essere piu' grande del picco raw (asintoto)
        %         p2 > 0, p3 in [0.5, 10]
        lb_rise = [min(Apeak_raw, 0)*2,  0.05,  0.5];
        ub_rise = [max(Apeak_raw, 0)*5,  50,    10 ];

        try
            pRise = lsqcurvefit(riseModel, p0_rise, t_rise, y_rise, ...
                                lb_rise, ub_rise, opts_lsq);
            Apeak_fit = riseModel(pRise, t_peak);  % ampiezza al picco dal fit
        catch
            pRise     = [NaN NaN NaN];
            Apeak_fit = Apeak_raw;
        end

        % ================================================================
        % FIT 2: discesa logistica con condizione al contorno
        %   y(t_decay=0) = Apeak_fit  ESATTAMENTE
        %
        %   Formula: y = Apeak_fit * (1 + exp(-k*t_half))
        %                           / (1 + exp(k*(t_decay - t_half)))
        %
        %   A t_decay=0:  y = Apeak_fit * (1+exp(-k*t_half))/(1+exp(-k*t_half))
        %                   = Apeak_fit  ✓
        %   Al flesso (t_decay=t_half): y ~ Apeak_fit/2  (come logistica classica)
        % ================================================================
        decayModel = @(p, td) Apeak_fit .* (1 + exp(-p(1).*p(2))) ...
                              ./ (1 + exp(p(1).*(td - p(2))));

        t_end_decay = t_decay(end);
        k_0    = 4 / max(t_end_decay, 1);     % pendenza iniziale
        th_0   = t_end_decay / 3;             % flesso a 1/3 della discesa

        p0_decay = [k_0, th_0];
        lb_decay = [0.01,  0            ];
        ub_decay = [50,    t_end_decay*2 ];

        try
            pDecay = lsqcurvefit(decayModel, p0_decay, t_decay, y_decay, ...
                                 lb_decay, ub_decay, opts_lsq);
        catch
            pDecay = [NaN NaN];
        end

        % ================================================================
        % Costruzione traccia fittata
        % ================================================================
        avatar = zeros(nTime, 1);

        if ~any(isnan(pRise)) && ~any(isnan(pDecay))
            avatar(riseMask)  = riseModel(pRise,  t_stim(riseMask));
            avatar(decayMask) = decayModel(pDecay, t_stim(decayMask) - t_peak);

            % Picco VERO della curva fittata (massimo della salita fittata)
            t_dense_rise  = linspace(0, t_peak, 1000)';
            y_dense_rise  = riseModel(pRise, t_dense_rise);
            [Apeak_plot, idx_plot] = max(y_dense_rise);
            t_peak_plot   = t_dense_rise(idx_plot);

            % R^2 calcolato sul post-stimolo (traccia raw)
            y_data = y_raw(postIdx);
            y_fit  = avatar(postIdx);
            SS_res = sum((y_data - y_fit).^2, 'omitnan');
            SS_tot = sum((y_data - mean(y_data,'omitnan')).^2, 'omitnan');
            R2 = 1 - SS_res / max(SS_tot, eps);
        else
            % Fallback: usa smooth
            avatar(postIdx) = y_smooth(postIdx);
            R2 = NaN;
            t_peak_plot = t_peak;
            Apeak_plot  = Apeak_raw;
            fprintf('  Stimolo %s, Animale %d: fit fallito, uso traccia smoothata.\n', ...
                    stimLabel, iAnim);
        end

        R2_stim(iAnim)     = R2;
        params_stim{iAnim} = [pRise, pDecay];

        % Salva risultati strutturati per questo animale
        animalResult.p1        = pRise(1);       % ampiezza asintotica salita
        animalResult.p2        = pRise(2);       % rate constant salita
        animalResult.p3        = pRise(3);       % shape salita
        animalResult.k         = pDecay(1);      % pendenza discesa logistica
        animalResult.t_half    = pDecay(2);      % flesso discesa (secondi da t_peak)
        animalResult.t_peak    = t_peak_plot;    % tempo al picco (dal fit, s)
        animalResult.A_peak    = Apeak_plot;     % ampiezza al picco (dal fit)
        animalResult.R2        = R2;
        animalResult.avatar    = avatar;         % traccia fittata [nTime x 1]
        animalResult.t         = t_stim;         % vettore tempo [nTime x 1]
        animalResult.y_raw     = y_raw;          % dati grezzi baseline-corretti
        stimResults(iAnim)     = animalResult;

        % ================================================================
        % Plot
        % ================================================================
        ax = subplot(nRows, nCols, iAnim, 'Parent', figH);
        hold(ax, 'on');

        % Dati grezzi
        plot(ax, t_stim, y_raw, 'Color', [0.75 0.75 0.75], 'LineWidth', 1);

        % Fit colorato (solo post-stimolo)
        t_fit_plot = t_stim(postIdx);
        plot(ax, t_fit_plot, avatar(postIdx), '-', ...
             'Color', animalColors(iAnim, :), 'LineWidth', 2);

        % Picco trovato (sul fit, non sui dati)
        Apeak_smooth = y_post_smooth(peak_idx_local);
        plot(ax, t_peak_plot, Apeak_plot, 'v', ...
             'MarkerFaceColor', animalColors(iAnim,:), ...
             'MarkerEdgeColor', 'k', 'MarkerSize', 6);

        xline(ax, 0, '--k', 'Alpha', 0.4);

        xlabel(ax, 'Tempo (s)');
        ylabel(ax, '\DeltaF/F o a.u.');
        title(ax, sprintf('Animale %d | R^2 = %.3f', iAnim, R2));
        box(ax, 'off');

        % Annotazione parametri
        if ~any(isnan([pRise, pDecay]))
            txt = sprintf('p1=%.2g  p2=%.2g  p3=%.2g\nk=%.2g  t_{1/2}=%.2g\nt_{pk}=%.3g  A_{pk}=%.3g', ...
                pRise(1), pRise(2), pRise(3), pDecay(1), pDecay(2), ...
                t_peak_plot, Apeak_plot);
            text(ax, 0.02, 0.98, txt, 'Units', 'normalized', ...
                'VerticalAlignment', 'top', 'FontSize', 7, ...
                'Interpreter', 'tex');
        end

        hold(ax, 'off');
    end

    sgtitle(figH, ['Stimolo: ' stimLabel], 'FontSize', 14, 'FontWeight', 'bold');

    allR2{iStim}     = R2_stim;
    allParams{iStim} = params_stim;

    % Assegna variabile con nome del stimolo nel workspace
    assignin('base', ['fitResults_' stimLabel], stimResults);
    fprintf('Variabile fitResults_%s creata (%d animali).\n', stimLabel, nAnimals);

    exportgraphics(figH, ['fit_' stimLabel '.pdf'], 'ContentType', 'vector');
    fprintf('Figura salvata: fit_%s.pdf\n', stimLabel);
end

%% ---- Figura riassuntiva R^2 ----
figSumm = figure('Name', 'R^2 per stimolo', 'NumberTitle', 'off', 'Color', 'w');
ax_s = axes(figSumm);
hold(ax_s, 'on');

jitter_amount = 0.12;

for iStim = 1:nStim
    R2_vals = allR2{iStim};
    nA      = numel(R2_vals);
    xPos    = iStim + (rand(1, nA) - 0.5) * jitter_amount * 2;

    scatter(ax_s, xPos, R2_vals, 70, ...
        'filled', 'MarkerFaceColor', animalColors(iStim,:), ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.8, ...
        'MarkerFaceAlpha', 0.8);

    % Mediana
    med = median(R2_vals, 'omitnan');
    plot(ax_s, [iStim-0.2, iStim+0.2], [med, med], '-k', 'LineWidth', 2);
end

set(ax_s, 'XTick', 1:nStim, 'XTickLabel', stimNames, 'FontSize', 12);
xlim(ax_s, [0.5, nStim + 0.5]);
ylim(ax_s, [0, 1.05]);
ylabel(ax_s, 'R^2');
xlabel(ax_s, 'Stimolo visivo');
title(ax_s, 'Bonta'' del fit logistico per stimolo e animale');
box(ax_s, 'off');
grid(ax_s, 'on');
grid(ax_s, 'minor');

exportgraphics(figSumm, 'R2_summary.pdf', 'ContentType', 'vector');
fprintf('Figura riassuntiva salvata: R2_summary.pdf\n');

%% fit_GRID5_expdecay.m
%
% Fit tripartito per stimolo GRID5, con decay ESPONENZIALE.
% Struttura analoga al codice di riferimento (fit 830nm).
%
% Modello:
%   t < 0            : M(t) = 0
%   0 <= t <= t_peak : M(t) = p1 * (1 - exp(-p2*t))^p3       [salita]
%   t >  t_peak      : M(t) = Apeak_fit * exp(-t_decay / tau) [decay esponenziale]
%
% Variabili richieste nel workspace:
%   t     - vettore tempo [nTime x 1], con t=0 = inizio stimolo
%   GRID5 - matrice [nTime x nAnimals]
%
% Output nel workspace:
%   fitResults_GRID5  - struct array [1 x nAnimals] con tutti i risultati
%
% Campi di fitResults_GRID5(i):
%   .p1, .p2, .p3   parametri salita
%   .tau            costante di tempo del decay esponenziale (s)
%   .Apeak_fit      ampiezza al picco ricavata dal fit di salita
%   .t_peak         tempo al picco (s), trovato su traccia smoothata
%   .t_peak_fit     tempo al picco della curva fittata (s)
%   .A_peak_fit     ampiezza della curva fittata al suo picco
%   .R2             R² calcolato sul post-stimolo (vs traccia raw)
%   .avatar         traccia fittata completa [nTime x 1]
%   .t              vettore tempo [nTime x 1]
%   .y_raw          dati grezzi baseline-corretti [nTime x 1]
%   .y_smooth       traccia smoothata baseline-corretta [nTime x 1]

%% ---- Configurazione ----

% Finestra di smoothing SIMMETRICA per trovare il picco (campioni prima e dopo)
smoothWin = 1;

% Finestra temporale in cui cercare il picco (secondi)
% Modifica questi valori in base alla tua finestra biologica attesa
tPeakSearch_start = 0;    % s dopo stimolo
tPeakSearch_end   = 10;   % s dopo stimolo

% Opzioni ottimizzatore (identiche al codice di riferimento)
opts_lsq = optimoptions('lsqcurvefit', ...
    'Display',                'off', ...
    'MaxFunctionEvaluations', 5000, ...
    'MaxIterations',          2000, ...
    'FunctionTolerance',      1e-10, ...
    'StepTolerance',          1e-10);

% Colori animali
nAnimals     = size(GRID5, 2);
nTime        = size(GRID5, 1);
animalColors = lines(nAnimals);

% Vettore tempo per GRID5
t_stim = t(1:nTime);

% Indici
baselineIdx  = t_stim < 0;
postIdx      = t_stim >= 0;
peakWinIdx   = (t_stim >= tPeakSearch_start) & (t_stim <= tPeakSearch_end);

%% ---- Storage ----
clear fitResults_GRID5   % evita conflitti con run precedenti

%% ---- Figura ----
nCols = ceil(sqrt(nAnimals));
nRows = ceil(nAnimals / nCols);
figH  = figure('Name', 'Fit GRID5 - Exp Decay', ...
               'NumberTitle', 'off', ...
               'Color', 'w', ...
               'Units', 'normalized', ...
               'OuterPosition', [0 0 1 1]);

%% ---- Loop animali ----
for iAnim = 1:nAnimals

    % --- Dati grezzi, correzione baseline ---
    y_raw   = GRID5(:, iAnim);
    y_raw   = y_raw - mean(y_raw(baselineIdx), 'omitnan');

    % --- Smooth SIMMETRICO per ricerca picco (no ritardo di fase) ---
    y_smooth = movmean(y_raw, [smoothWin smoothWin], 'omitnan');
    y_smooth = y_smooth - mean(y_smooth(baselineIdx), 'omitnan');

    % --- Ricerca picco sulla traccia smoothata nella finestra peakWin ---
    % (come nel codice di riferimento con peakWin_830)
    tSearch       = t_stim(peakWinIdx);
    [~, idxLocal] = max(y_smooth(peakWinIdx));
    t_peak        = tSearch(idxLocal);

    % Indice globale del picco
    [~, peakIdx_global] = min(abs(t_stim - t_peak));
    Apeak_smooth        = y_smooth(peakIdx_global);   % valore smooth al picco
    Apeak_raw           = y_raw(peakIdx_global);      % valore raw al picco

    % --- Maschere salita e discesa ---
    riseMask  = (t_stim >= 0) & (t_stim <= t_peak);
    decayMask =  t_stim >= t_peak;

    t_rise  = t_stim(riseMask);
    y_rise  = y_raw(riseMask);
    t_decay = t_stim(decayMask) - t_peak;   % riscalato da 0 (come ref)
    y_decay = y_raw(decayMask);

    % ================================================================
    % FIT 1: salita  ->  p1*(1-exp(-p2*t))^p3
    % (identico al codice di riferimento)
    % ================================================================
    riseModel = @(p, t) p(1) .* (1 - exp(-p(2) .* t)) .^ p(3);

    p0_rise = [Apeak_raw,          1.5,  2  ];
    lb_rise = [min(5*Apeak_raw,0), 0.25, 1  ];
    ub_rise = [max(5*Apeak_raw,0), 20,   10 ];

    try
        pRise     = lsqcurvefit(riseModel, p0_rise, t_rise, y_rise, ...
                                lb_rise, ub_rise, opts_lsq);
        Apeak_fit = riseModel(pRise, t_rise(end));   % ampiezza al picco dal fit
    catch
        pRise     = [NaN NaN NaN];
        Apeak_fit = Apeak_raw;
    end

    % ================================================================
    % FIT 2: decay esponenziale
    %   y = Apeak_fit * exp(-t_decay / tau)
    %   Continuita' di valore garantita: a t_decay=0 -> Apeak_fit ✓
    % ================================================================
    decayModel = @(tau, t) Apeak_fit .* exp(-t ./ tau);

    tau0    = 5;    % stima iniziale (s) — regola se necessario
    lb_tau  = 0.5;
    ub_tau  = 60;

    try
        tauFit = lsqcurvefit(decayModel, tau0, t_decay, y_decay, ...
                             lb_tau, ub_tau, opts_lsq);
    catch
        tauFit = NaN;
    end

    % ================================================================
    % Costruzione traccia fittata (avatar)
    % ================================================================
    avatar = zeros(nTime, 1);

    if ~any(isnan(pRise)) && ~isnan(tauFit)

        avatar(riseMask)  = riseModel(pRise,  t_stim(riseMask));
        avatar(decayMask) = decayModel(tauFit, t_stim(decayMask) - t_peak);

        % Picco VERO della curva fittata (su griglia densa)
        t_dense      = linspace(0, t_peak, 2000)';
        y_dense      = riseModel(pRise, t_dense);
        [A_peak_fit, idx_pk] = max(y_dense);
        t_peak_fit   = t_dense(idx_pk);

        % R² sul post-stimolo vs dati RAW (come codice di riferimento)
        y_data = y_raw(postIdx);
        y_fit  = avatar(postIdx);
        SS_res = sum((y_data - y_fit).^2, 'omitnan');
        SS_tot = sum((y_data - mean(y_data,'omitnan')).^2, 'omitnan');
        if SS_tot > 0
            R2 = 1 - SS_res / SS_tot;
        else
            R2 = NaN;
        end

    else
        % Fallback: smooth
        avatar(postIdx) = y_smooth(postIdx);
        R2          = NaN;
        t_peak_fit  = t_peak;
        A_peak_fit  = Apeak_raw;
        fprintf('  GRID5, Animale %d: fit fallito, uso traccia smoothata.\n', iAnim);
    end

    % ================================================================
    % Salva risultati strutturati
    % ================================================================
    fitResults_GRID5(iAnim).p1          = pRise(1);     % ampiezza asintotica salita
    fitResults_GRID5(iAnim).p2          = pRise(2);     % rate constant salita
    fitResults_GRID5(iAnim).p3          = pRise(3);     % shape salita
    fitResults_GRID5(iAnim).tau         = tauFit;       % costante decay (s)
    fitResults_GRID5(iAnim).Apeak_fit   = Apeak_fit;    % ampiezza al picco (fit salita)
    fitResults_GRID5(iAnim).t_peak      = t_peak;       % t picco da smooth (s)
    fitResults_GRID5(iAnim).t_peak_fit  = t_peak_fit;   % t picco della curva fittata (s)
    fitResults_GRID5(iAnim).A_peak_fit  = A_peak_fit;   % ampiezza picco curva fittata
    fitResults_GRID5(iAnim).R2          = R2;
    fitResults_GRID5(iAnim).avatar      = avatar;       % [nTime x 1]
    fitResults_GRID5(iAnim).t           = t_stim;       % [nTime x 1]
    fitResults_GRID5(iAnim).y_raw       = y_raw;        % [nTime x 1]
    fitResults_GRID5(iAnim).y_smooth    = y_smooth;     % [nTime x 1]

    % ================================================================
    % Plot
    % ================================================================
    ax = subplot(nRows, nCols, iAnim, 'Parent', figH);
    hold(ax, 'on');

    % Dati grezzi
    plot(ax, t_stim, y_raw, 'Color', [0.75 0.75 0.75], 'LineWidth', 1);

    % Fit (solo post-stimolo)
    plot(ax, t_stim(postIdx), avatar(postIdx), '-', ...
         'Color', animalColors(iAnim,:), 'LineWidth', 2);

    % Picco della curva fittata
    plot(ax, t_peak_fit, A_peak_fit, 'v', ...
         'MarkerFaceColor', animalColors(iAnim,:), ...
         'MarkerEdgeColor', 'k', 'MarkerSize', 7);

    xline(ax, 0, '--k', 'Alpha', 0.4);

    xlabel(ax, 'Tempo (s)');
    ylabel(ax, '\DeltaF/F o a.u.');
    title(ax, sprintf('Animale %d | R^2 = %.3f', iAnim, R2));
    box(ax, 'off');

    % Annotazione parametri
    if ~any(isnan(pRise)) && ~isnan(tauFit)
        txt = sprintf('p1=%.2g  p2=%.2g  p3=%.2g\n\\tau=%.2g s\nt_{pk}=%.2g s  A_{pk}=%.2g', ...
            pRise(1), pRise(2), pRise(3), tauFit, t_peak_fit, A_peak_fit);
        text(ax, 0.02, 0.98, txt, ...
            'Units', 'normalized', ...
            'VerticalAlignment', 'top', ...
            'FontSize', 8, ...
            'Interpreter', 'tex');
    end

    hold(ax, 'off');
end

sgtitle(figH, 'GRID5 — Fit tripartito con decay esponenziale', ...
        'FontSize', 14, 'FontWeight', 'bold');

%exportgraphics(figH, 'fit_GRID5_expdecay.pdf', 'ContentType', 'vector');
fprintf('Figura salvata: fit_GRID5_expdecay.pdf\n');

%% ---- Riepilogo a schermo ----
fprintf('\n===== Riepilogo fit GRID5 =====\n');
fprintf('%-8s  %6s  %6s  %6s  %8s  %8s  %6s\n', ...
        'Animale','p1','p2','p3','tau (s)','t_pk (s)','R²');
for iAnim = 1:nAnimals
    r = fitResults_GRID5(iAnim);
    fprintf('%-8d  %6.2f  %6.2f  %6.2f  %8.2f  %8.3f  %6.3f\n', ...
            iAnim, r.p1, r.p2, r.p3, r.tau, r.t_peak_fit, r.R2);
end
fprintf('================================\n');

%% ---- Come usare i risultati ----
%
% Parametri numerici per il topo i:
%   fitResults_GRID5(i).tau        -> costante di decay (s)
%   fitResults_GRID5(i).t_peak_fit -> tempo al picco (s)
%   fitResults_GRID5(i).A_peak_fit -> ampiezza al picco
%   fitResults_GRID5(i).R2         -> bonta' del fit
%   fitResults_GRID5(i).p1/p2/p3   -> parametri salita
%
% Traccia fittata (avatar) del topo i, da plottare altrove:
%   plot(fitResults_GRID5(i).t, fitResults_GRID5(i).avatar)
%
% Estrarre tutti i tau in un vettore:
%   tau_all  = [fitResults_GRID5.tau];
%
% Estrarre tutti gli R² in un vettore:
%   R2_all   = [fitResults_GRID5.R2];
%
% Estrarre tutti i t_peak in un vettore:
%   tpk_all  = [fitResults_GRID5.t_peak_fit];