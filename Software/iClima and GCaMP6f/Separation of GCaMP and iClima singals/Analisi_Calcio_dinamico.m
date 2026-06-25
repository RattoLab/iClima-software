% This script was used to extract GCaMP6f responses from double positive
% cells. Moreover it identified 'unmasked' GCaMP6f responses
% =========================================================================
% ClCa_session2_stimAndSmooth.m
% =========================================================================
%
% DESCRIZIONE:
%   Seconda sessione di analisi sulla struttura ClCa.
%   Questo script esegue due operazioni principali:
%
%   1) ASSEGNAZIONE DEL TIPO DI STIMOLO
%      Per ogni riga di ClCa, analizza il nome del file (FileName_830)
%      per identificare il tipo di stimolo somministrato. I tipi riconosciuti
%      sono:
%        - GRID con orientamento specifico (0, 45, 90, 135, 180, 270, 315)
%        - CYCLEGRID (stimolo ciclico a due orientamenti, es. 0-90)
%      Se lo stimolo non viene riconosciuto, la riga viene segnalata in
%      console e il campo Stimulus rimane vuoto.
%
%   2) CALCOLO DEL FRAME PERIOD E SMOOTHING CAUSALE
%      Per ogni riga di ClCa, a 830 nm e a 920 nm:
%        - Calcola il frame period (mediana di diff(Time))
%        - Esegue uno smooth causale (moving average sul passato) di
%          G_sub_dark e R_sub_dark, usando una finestra temporale fissa
%          definita dall'operatore in secondi. La finestra viene convertita
%          in numero di frames in base al frame period di ciascuna
%          acquisizione, garantendo che acquisizioni a frequenze diverse
%          vengano trattate sulla stessa scala temporale.
%
% PARAMETRI EDITABILI (Sezione 0):
%   smooth_window_sec : finestra temporale per lo smooth [secondi]
%
% INPUT:
%   La struttura ClCa deve essere presente nel workspace (caricata
%   dallo script build_ClCa.m o da un file .mat salvato in precedenza).
%
% OUTPUT:
%   Aggiorna ClCa nel workspace aggiungendo i campi:
%     - Stimulus           : tipo di stimolo assegnato (stringa)
%     - FramePeriod_830/920: frame period in secondi
%     - G_smooth_830/920   : canale verde smoothato (causale)
%     - R_smooth_830/920   : canale rosso smoothato (causale)
%
% =========================================================================

fprintf('=========================================================\n');
fprintf('   ClCa SESSION 2 - Stimoli e Smoothing causale          \n');
fprintf('=========================================================\n\n');

% Verifica che ClCa sia presente nel workspace
if ~exist('ClCa','var') || isempty(ClCa)
    % Se non è nel workspace, chiede di caricarla
    [fname_clca, fdir_clca] = uigetfile('*.mat', ...
        'ClCa non trovata nel workspace. Seleziona il file .mat');
    if isequal(fname_clca, 0)
        error('Nessun file selezionato. Script terminato.');
    end
    tmp = load(fullfile(fdir_clca, fname_clca), 'ClCa');
    if ~isfield(tmp, 'ClCa')
        error('Il file selezionato non contiene una variabile "ClCa".');
    end
    ClCa = tmp.ClCa;
    fprintf('ClCa caricata dal file: %s\n\n', fname_clca);
end

fprintf('ClCa contiene %i righe.\n\n', numel(ClCa));


% =========================================================================
%% SEZIONE 0 — PARAMETRI EDITABILI
% =========================================================================
% Finestra temporale per lo smooth causale, in secondi.
% Verrà convertita in numero di frames per ogni acquisizione
% in base al suo frame period, così acquisizioni a diverse frequenze
% vengono smoothate sulla stessa scala temporale.
%
% Esempio: 0.5 secondi equivale a:
%   - 1 frame  se acquisito a 2  Hz  (frame period 0.5 s)
%   - 15 frame se acquisito a 30 Hz  (frame period 0.033 s)

smooth_window_sec = 1;   % <-- MODIFICA QUI la finestra in secondi

fprintf('Finestra di smooth: %.3f secondi\n\n', smooth_window_sec);


% =========================================================================
%% SEZIONE 1 — ASSEGNAZIONE DEL TIPO DI STIMOLO
% =========================================================================
% Per ogni riga di ClCa, il tipo di stimolo viene determinato analizzando
% il nome del file a 830 nm (FileName_830). Stesso stimolo per entrambe
% le lunghezze d'onda (sono acquisizioni dello stesso replicato biologico).
%
% Gerarchia di controllo:
%   1. Prima controlla CYCLEGRID (ha la precedenza perché contiene anche
%      numeri di griglia nel nome)
%   2. Poi controlla GRID con orientamento specifico

fprintf('--- Sezione 1: Assegnazione stimoli ---\n');

nUnmatched = 0;   % contatore righe senza stimolo riconosciuto

for i = 1 : numel(ClCa)

    fname = ClCa(i).FileName_830;   % usa il nome del file 830 come riferimento

    % Tenta l'assegnazione dello stimolo
    [stimLabel, stimFound] = assignStimulus(fname);

    if stimFound
        ClCa(i).Stimulus = stimLabel;
        fprintf('  Riga %3i | %-30s -> Stimulus: %s\n', i, fname, stimLabel);
    else
        % Stimolo non riconosciuto: segnala in console, lascia campo vuoto
        ClCa(i).Stimulus = '';
        nUnmatched = nUnmatched + 1;
        fprintf('  [NON TROVATO] Riga %3i | Mouse: %s | File: %s\n', ...
            i, ClCa(i).Mouse, fname);
    end

end

fprintf('\nAssegnazione stimoli completata.\n');
fprintf('  Righe con stimolo assegnato : %i / %i\n', numel(ClCa) - nUnmatched, numel(ClCa));
fprintf('  Righe senza corrispondenza  : %i\n\n', nUnmatched);


% =========================================================================
%% SEZIONE 2 — FRAME PERIOD E SMOOTHING CAUSALE
% =========================================================================
% Per ogni riga di ClCa, e per ciascuna delle due lunghezze d'onda (830
% e 920), lo script:
%   a) Calcola il frame period come mediana di diff(Time).
%      Viene usata la mediana (robusta agli outlier) anziché la media.
%   b) Converte la finestra temporale smooth_window_sec in numero di
%      frames: win_frames = round(smooth_window_sec / frame_period)
%      con un minimo di 1 frame.
%   c) Applica uno smooth causale (moving average) a G_sub_dark e
%      R_sub_dark tramite la funzione causalSmooth().
%      "Causale" significa che la media è calcolata solo sui frame
%      precedenti (incluso quello corrente), mai sui frame futuri.
%      Questo evita che la risposta sembri anticipare lo stimolo.

fprintf('--- Sezione 2: Frame period e smoothing causale ---\n');

for i = 1 : numel(ClCa)

    fprintf('  Elaboro riga %i/%i (Mouse: %s, Field: %s)...\n', ...
        i, numel(ClCa), ClCa(i).Mouse, ClCa(i).Field);

    % ---- 830 nm ----------------------------------------------------------
    if ~isempty(ClCa(i).Time_830)

        % a) Frame period a 830 nm
        fp830 = median(diff(ClCa(i).Time_830));
        ClCa(i).FramePeriod_830 = fp830;

        % b) Numero di frames corrispondente alla finestra temporale
        win_frames_830 = max(1, round(smooth_window_sec / fp830));
        fprintf('     830nm: fp=%.4fs | finestra smooth=%i frames\n', ...
            fp830, win_frames_830);

        % c) Smooth causale di G e R a 830 nm
        if ~isempty(ClCa(i).G_sub_dark_830)
            ClCa(i).G_smooth_830 = causalSmooth(ClCa(i).G_sub_dark_830, win_frames_830);
        else
            ClCa(i).G_smooth_830 = [];
            warning('G_sub_dark_830 vuoto per riga %i', i);
        end

        if ~isempty(ClCa(i).R_sub_dark_830)
            ClCa(i).R_smooth_830 = causalSmooth(ClCa(i).R_sub_dark_830, win_frames_830);
        else
            ClCa(i).R_smooth_830 = [];
            warning('R_sub_dark_830 vuoto per riga %i', i);
        end

    else
        % Time_830 assente: impossibile calcolare frame period
        ClCa(i).FramePeriod_830 = NaN;
        ClCa(i).G_smooth_830    = [];
        ClCa(i).R_smooth_830    = [];
        warning('Time_830 vuoto per riga %i: smooth non eseguito.', i);
    end

    % ---- 920 nm ----------------------------------------------------------
    if ~isempty(ClCa(i).Time_920)

        % a) Frame period a 920 nm
        fp920 = median(diff(ClCa(i).Time_920));
        ClCa(i).FramePeriod_920 = fp920;

        % b) Numero di frames corrispondente alla finestra temporale
        win_frames_920 = max(1, round(smooth_window_sec / fp920));
        fprintf('     920nm: fp=%.4fs | finestra smooth=%i frames\n', ...
            fp920, win_frames_920);

        % c) Smooth causale di G e R a 920 nm
        if ~isempty(ClCa(i).G_sub_dark_920)
            ClCa(i).G_smooth_920 = causalSmooth(ClCa(i).G_sub_dark_920, win_frames_920);
        else
            ClCa(i).G_smooth_920 = [];
            warning('G_sub_dark_920 vuoto per riga %i', i);
        end

        if ~isempty(ClCa(i).R_sub_dark_920)
            ClCa(i).R_smooth_920 = causalSmooth(ClCa(i).R_sub_dark_920, win_frames_920);
        else
            ClCa(i).R_smooth_920 = [];
            warning('R_sub_dark_920 vuoto per riga %i', i);
        end

    else
        ClCa(i).FramePeriod_920 = NaN;
        ClCa(i).G_smooth_920    = [];
        ClCa(i).R_smooth_920    = [];
        warning('Time_920 vuoto per riga %i: smooth non eseguito.', i);
    end

end

fprintf('\nSmoothing completato.\n\n');


% =========================================================================
%% SEZIONE 3 — SALVATAGGIO
% =========================================================================

risposta_salva = questdlg('Vuoi salvare la struttura ClCa aggiornata?', ...
    'Salvataggio', 'Sì', 'No (mantieni solo in workspace)', 'Sì');

if strcmp(risposta_salva, 'Sì')
    defaultName = ['ClCa_session2_', datestr(now, 'yyyymmdd_HHMM'), '.mat']; %#ok<TNOW1,DATST>
    [saveName, saveDir] = uiputfile('*.mat', ...
        'Salva ClCa aggiornata', defaultName);
    if ~isequal(saveName, 0)
        save(fullfile(saveDir, saveName), 'ClCa');
        fprintf('ClCa salvata in:\n  %s\n', fullfile(saveDir, saveName));
    else
        fprintf('Salvataggio annullato. ClCa disponibile nel workspace.\n');
    end
else
    fprintf('ClCa disponibile nel workspace (non salvata su file).\n');
end

fprintf('\n=========================================================\n');
fprintf('               Sessione 2 completata.                    \n');
fprintf('=========================================================\n');


% =========================================================================


%% =========================================================================
% ClCa_session3_avatarAndScaling.m  (G/R version)
% =========================================================================

fprintf('=========================================================\n');
fprintf('   ClCa SESSION 3 - Avatar e Scaling 830->920  (G/R)     \n');
fprintf('=========================================================\n\n');

if ~exist('ClCa','var') || isempty(ClCa)
    [fname_clca, fdir_clca] = uigetfile('*.mat', ...
        'ClCa non trovata nel workspace. Seleziona il file .mat');
    if isequal(fname_clca, 0)
        error('Nessun file selezionato. Script terminato.');
    end
    tmp = load(fullfile(fdir_clca, fname_clca), 'ClCa');
    if ~isfield(tmp, 'ClCa')
        error('Il file selezionato non contiene una variabile "ClCa".');
    end
    ClCa = tmp.ClCa;
    fprintf('ClCa caricata: %i righe.\n\n', numel(ClCa));
end


% =========================================================================
%% SEZIONE 0 — PARAMETRI EDITABILI
% =========================================================================

stimToAnalyze = {'Grid_0','Grid_45','Grid_90','Grid_180','Grid_315'};   % <-- MODIFICA QUI

peak_search_end_sec = 7;        % <-- MODIFICA QUI  finestra ricerca [s]
scaling_window_end  = 4.5;      % <-- MODIFICA QUI  finestra scaling [s]

% --- Soglie classificazione ROI ------------------------------------------
%
% alpha_threshold : ROI con scaleFactor < soglia → AbnormalFit (asse rosso)
%   se non già classificata come solo-calcio o grigia.
alpha_threshold = 0.15;         % <-- MODIFICA QUI

% calcium_only_threshold : soglia sul picco POSITIVO di dG/R a 920nm.
%   Espresso in frazione normalizzata (0.20 = 20% sopra baseline).
calcium_only_threshold = 0.20;  % <-- MODIFICA QUI

% min_negative_peak_830 : ampiezza minima della valle a 830nm (valore
%   assoluto). Se abs(minimo dG/R a 830nm) < questa soglia, iClima non
%   risponde significativamente → condizione necessaria per solo-calcio.
%   Espresso in frazione normalizzata (0.02 = 2%).
min_negative_peak_830 = 0.01;   % <-- MODIFICA QUI

fprintf('Parametri classificazione:\n');
fprintf('  alpha_threshold        = %.3f\n', alpha_threshold);
fprintf('  calcium_only_threshold = %.3f\n', calcium_only_threshold);
fprintf('  min_negative_peak_830  = %.3f\n\n', min_negative_peak_830);


% =========================================================================
%% SEZIONE 1 — SELEZIONE STIMOLO E METODO DI SCALING
% =========================================================================

if ischar(stimToAnalyze)
    stimToAnalyze = {stimToAnalyze};
end

fprintf('Stimoli selezionati: %s\n', strjoin(stimToAnalyze, ', '));
rowsWithStim = find(ismember({ClCa.Stimulus}, stimToAnalyze));
if isempty(rowsWithStim)
    error('Nessuna riga di ClCa ha gli stimoli specificati.');
end
fprintf('Righe trovate: %i\n\n', numel(rowsWithStim));

scalingMethod = questdlg( ...
    'Quale metodo di scaling vuoi usare?', 'Metodo di scaling', ...
    'Peak-based', 'Residual minimization', 'Peak-based');
if isempty(scalingMethod)
    error('Nessun metodo selezionato. Script terminato.');
end
fprintf('Metodo di scaling scelto: %s\n\n', scalingMethod);

opts_lsq = optimoptions('lsqcurvefit', 'Display', 'off');


% =========================================================================
%% SEZIONE 2 — LOOP SULLE RIGHE VALIDE
% =========================================================================

fprintf('--- Inizio elaborazione righe ---\n\n');c

for iRow = rowsWithStim

    fprintf('Elaboro riga %i | Mouse: %s | Field: %s | File830: %s\n', ...
        iRow, ClCa(iRow).Mouse, ClCa(iRow).Field, ClCa(iRow).FileName_830);

    nROI_830 = size(ClCa(iRow).G_smooth_830, 3);
    nROI_920 = size(ClCa(iRow).G_smooth_920, 3);
    if nROI_830 ~= nROI_920
        fprintf('  [SALTATA] ROI non coerenti: 830nm=%i, 920nm=%i\n\n', nROI_830, nROI_920);
        continue
    end
    nROIs = nROI_830;
    fprintf('  ROI totali: %i\n', nROIs);

    Ratio_830 = ClCa(iRow).G_smooth_830 ./ ClCa(iRow).R_smooth_830;
    Ratio_920 = ClCa(iRow).G_smooth_920 ./ ClCa(iRow).R_smooth_920;

    Time_830 = ClCa(iRow).Time_830(:);
    Time_920 = ClCa(iRow).Time_920(:);

    baselineIdx_830 = Time_830 < 0;
    baselineIdx_920 = Time_920 < 0;
    peakWin_830 = (Time_830 >= 0) & (Time_830 <= peak_search_end_sec);
    peakWin_920 = (Time_920 >= 0) & (Time_920 <= peak_search_end_sec);

    nTime_830 = length(Time_830);
    nTime_920 = length(Time_920);

    % Pre-allocazione
    avatar830_all      = nan(nTime_830, nROIs);
    mean830_all        = nan(nTime_830, nROIs);
    mean920_all        = nan(nTime_920, nROIs);
    scaledFit_all      = nan(nTime_920, nROIs);
    residual_all       = nan(nTime_920, nROIs);
    scaleFactor_all    = nan(nROIs, 1);
    tau830_all         = nan(nROIs, 1);
    fitParams_all      = nan(nROIs, 4);
    r2_830             = nan(nROIs, 1);

    % Metriche per classificazione
    maxPositivePeak_920 = nan(nROIs, 1);  % picco positivo dG/R a 920nm
    minNegativePeak_830 = nan(nROIs, 1);  % valore assoluto della valle a 830nm

    for roi = 1:nROIs

        % == 830 nm ========================================================

        ratio830_roi = squeeze(Ratio_830(:, :, roi));
        ratio0_830   = mean(ratio830_roi(baselineIdx_830, :), 1);
        dRatio_830   = (ratio830_roi - ratio0_830) ./ ratio0_830;

        meanRaw_830    = mean(dRatio_830, 2, 'omitnan');
        meanRaw_830    = meanRaw_830 - mean(meanRaw_830(baselineIdx_830), 'omitnan');
        meanSmooth_830 = movmean(meanRaw_830, [15 0], 'omitnan');
        meanSmooth_830 = meanSmooth_830 - mean(meanSmooth_830(baselineIdx_830), 'omitnan');

        mean830_all(:, roi) = meanRaw_830;

        % Ampiezza della valle a 830nm (valore assoluto del minimo)
        % → misura quanto iClima risponde allo stimolo
        minNegativePeak_830(roi) = abs(min(meanSmooth_830(peakWin_830)));

        % Ricerca del MINIMO (valle iClima) per il fit
        tSearch_830      = Time_830(peakWin_830);
        [~, idxLocal]    = min(meanSmooth_830(peakWin_830));
        tPeak_830        = tSearch_830(idxLocal);
        [~, peakIdx_830] = min(abs(Time_830 - tPeak_830));
        Apeak_830        = meanSmooth_830(peakIdx_830);

        riseMask_830  = (Time_830 >= 0) & (Time_830 <= tPeak_830);
        decayMask_830 = (Time_830 >= tPeak_830);

        tRise_830  = Time_830(riseMask_830);
        yRise_830  = meanRaw_830(riseMask_830);
        tDecay_830 = Time_830(decayMask_830) - tPeak_830;
        yDecay_830 = meanRaw_830(decayMask_830);

        riseModel = @(p, t) p(1) .* (1 - exp(-p(2) .* t)) .^ p(3);
        p0_rise   = [Apeak_830, 1.5, 2];
        lb_rise   = [min(5*Apeak_830, 0), 0.25, 1];
        ub_rise   = [max(5*Apeak_830, 0), 20,   10];

        try
            pRise_830 = lsqcurvefit(riseModel, p0_rise, tRise_830, yRise_830, ...
                                    lb_rise, ub_rise, opts_lsq);
            Apeak_fit = riseModel(pRise_830, tRise_830(end));
        catch
            pRise_830 = [NaN NaN NaN];
            Apeak_fit = Apeak_830;
        end

        decayModel = @(tau, t) Apeak_fit .* exp(-t ./ tau);
        try
            tauFit_830 = lsqcurvefit(decayModel, 10, tDecay_830, yDecay_830, ...
                                     0.5, 7, opts_lsq);
        catch
            tauFit_830 = NaN;
        end

        avatar830 = zeros(nTime_830, 1);
        if ~any(isnan(pRise_830)) && ~isnan(tauFit_830)
            avatar830(riseMask_830)  = riseModel(pRise_830, Time_830(riseMask_830));
            avatar830(decayMask_830) = decayModel(tauFit_830, Time_830(decayMask_830) - tPeak_830);
        else
            avatar830(riseMask_830)  = meanSmooth_830(riseMask_830);
            avatar830(decayMask_830) = meanSmooth_830(decayMask_830);
            fprintf('  ROI %i: fit fallito, uso traccia smoothata.\n', roi);
        end

        avatar830_all(:, roi) = avatar830;
        tau830_all(roi)       = tauFit_830;
        fitParams_all(roi, :) = [pRise_830, tauFit_830];

        % R² del fit sulla traccia smoothata post-stimolo
        postStimMask_830 = Time_830 >= 0;
        y_data = meanSmooth_830(postStimMask_830);
        y_fit  = avatar830(postStimMask_830);
        SS_res = sum((y_data - y_fit).^2, 'omitnan');
        SS_tot = sum((y_data - mean(y_data, 'omitnan')).^2, 'omitnan');
        if SS_tot > 0
            r2_830(roi) = 1 - SS_res / SS_tot;
        end

        % == 920 nm ========================================================

        ratio920_roi = squeeze(Ratio_920(:, :, roi));
        ratio0_920   = mean(ratio920_roi(baselineIdx_920, :), 1);
        dRatio_920   = (ratio920_roi - ratio0_920) ./ ratio0_920;

        meanRaw_920    = mean(dRatio_920, 2, 'omitnan');
        meanSmooth_920 = movmean(meanRaw_920, [5 0], 'omitnan');
        meanSmooth_920 = meanSmooth_920 - mean(meanSmooth_920(baselineIdx_920), 'omitnan');

        mean920_all(:, roi) = meanSmooth_920;

        % Picco positivo a 920nm → proxy risposta calcio
        maxPositivePeak_920(roi) = max(meanSmooth_920(peakWin_920));

        % == SCALING =======================================================

        avatar830_on920 = interp1(Time_830, avatar830, Time_920, 'linear', NaN);
        y920       = meanSmooth_920;
        f830_on920 = avatar830_on920;

        switch scalingMethod
            case 'Peak-based'
                tSearch_920 = Time_920(peakWin_920);
                [~, idx920] = min(y920(peakWin_920));
                tPeak920    = tSearch_920(idx920);
                [~, peakIdx920] = min(abs(Time_920 - tPeak920));
                A920 = y920(peakIdx920);

                [~, idx830on920] = min(f830_on920(peakWin_920));
                tPeak830on920    = tSearch_920(idx830on920);
                [~, peakIdx830on920] = min(abs(Time_920 - tPeak830on920));
                A830on920 = f830_on920(peakIdx830on920);

                if A830on920 ~= 0
                    scaleFactor = A920 / A830on920;
                else
                    scaleFactor = NaN;
                end

            case 'Residual minimization'
                scalingWin_920 = (Time_920 >= 0) & (Time_920 <= scaling_window_end);
                y_fit_sc  = y920(scalingWin_920);
                f_fit_sc  = f830_on920(scalingWin_920);
                validMask = ~isnan(y_fit_sc) & ~isnan(f_fit_sc);
                y_fit_sc  = y_fit_sc(validMask);
                f_fit_sc  = f_fit_sc(validMask);
                if ~isempty(f_fit_sc) && any(f_fit_sc ~= 0)
                    costFun     = @(a) sum((y_fit_sc - a .* f_fit_sc).^2);
                    scaleFactor = fminsearch(costFun, 1);
                else
                    scaleFactor = NaN;
                end
        end

        scaleFactor_all(roi) = scaleFactor;

        if ~isnan(scaleFactor)
            scaledFit = scaleFactor .* f830_on920;
            residual  = y920 - scaledFit;
        else
            scaledFit = nan(nTime_920, 1);
            residual  = nan(nTime_920, 1);
        end

        scaledFit_all(:, roi) = scaledFit;
        residual_all(:, roi)  = residual;

        fprintf('  ROI %2i | alpha=%.3f | valle830=%.4f | picco920=%.4f\n', ...
            roi, scaleFactor, minNegativePeak_830(roi), maxPositivePeak_920(roi));

    end % fine loop ROI

    % ---- Classificazione ROI (dopo il loop, con tutte le metriche) ------
    %
    % GERARCHIA (in ordine di priorità):
    %
    %   1. SOLO CALCIO (blu):
    %      picco positivo a 920nm > calcium_only_threshold
    %      AND valle a 830nm < min_negative_peak_830
    %      → GCaMP risponde, iClima non risponde
    %
    %   2. ABNORMAL FIT (rosso):
    %      alpha < alpha_threshold (o NaN)
    %      AND non già classificata come solo-calcio
    %      → fit geometricamente sbagliato o scala impossibile
    %
    %   3. NORMALE (nero):
    %      tutto il resto con valle 830nm >= min_negative_peak_830
    %      → risposta iClima presente e fit accettabile
    %
    %   4. GRIGIA:
    %      nessuna delle categorie precedenti
    %      → segnale troppo debole o ambiguo per essere classificato

    isCalciumOnly = false(nROIs, 1);
    isAbnormalFit = false(nROIs, 1);
    isNormal      = false(nROIs, 1);
    % Le grigie sono quelle che rimangono false in tutti e tre i flag

    for roi = 1:nROIs
        hasPosCalcium  = maxPositivePeak_920(roi) > calcium_only_threshold;
        hasWeakiClima  = minNegativePeak_830(roi) < min_negative_peak_830;
        hasGoodiClima  = minNegativePeak_830(roi) >= min_negative_peak_830;
        hasBadAlpha    = isnan(scaleFactor_all(roi)) || scaleFactor_all(roi) < alpha_threshold;

        if hasPosCalcium %&& hasWeakiClima
            % iClima assente + calcio presente → solo calcio
            isCalciumOnly(roi) = true;
        elseif hasBadAlpha
            % fit abnormal (non già solo-calcio)
            isAbnormalFit(roi) = true;
        elseif hasGoodiClima
            % risposta iClima presente e fit ok → normale
            isNormal(roi) = true;
        end
        % se nessuno: rimane grigia (tutti false)
    end

    isGray = ~isCalciumOnly & ~isAbnormalFit & ~isNormal;

    fprintf('\n  Classificazione ROI:\n');
    fprintf('    Normali (nero)    : %i\n', sum(isNormal));
    fprintf('    Solo calcio (blu) : %i\n', sum(isCalciumOnly));
    fprintf('    Abnormal (rosso)  : %i\n', sum(isAbnormalFit));
    fprintf('    Grigie            : %i\n\n', sum(isGray));

    % ---- Salvataggio risultati ------------------------------------------
    fieldSuffix = strrep(ClCa(iRow).Stimulus, '_', '');
    res = ['Results_' fieldSuffix];

    ClCa(iRow).(res).Time_830             = Time_830;
    ClCa(iRow).(res).Time_920             = Time_920;
    ClCa(iRow).(res).dRatio_830           = mean830_all;
    ClCa(iRow).(res).dRatio_920           = mean920_all;
    ClCa(iRow).(res).Avatar_830           = avatar830_all;
    ClCa(iRow).(res).ScaledFit            = scaledFit_all;
    ClCa(iRow).(res).Residual             = residual_all;
    ClCa(iRow).(res).ScaleFactor          = scaleFactor_all;
    ClCa(iRow).(res).FitParams            = fitParams_all;
    ClCa(iRow).(res).Tau830               = tau830_all;
    ClCa(iRow).(res).ScalingMethod        = scalingMethod;
    ClCa(iRow).(res).R2_fit_830           = r2_830;
    ClCa(iRow).(res).IsNormal             = isNormal;
    ClCa(iRow).(res).IsCalciumOnly        = isCalciumOnly;
    ClCa(iRow).(res).IsAbnormalFit        = isAbnormalFit;
    ClCa(iRow).(res).IsGray               = isGray;
    % Retrocompatibilità con session4
    ClCa(iRow).(res).AbnormalAvatar       = isCalciumOnly | isAbnormalFit | isGray;
    ClCa(iRow).(res).AlphaThreshold       = alpha_threshold;
    ClCa(iRow).(res).CalciumOnlyThreshold = calcium_only_threshold;
    ClCa(iRow).(res).MinNegPeak830        = min_negative_peak_830;

    fprintf('  Risultati salvati in ClCa(%i).%s\n\n', iRow, res);

    % ---- Visualizzazione avatar ------------------------------------------
    % Ordine: normali → grigie → solo-calcio → abnormal
    norm_idx   = find(isNormal);
    gray_idx   = find(isGray);
    ca_idx     = find(isCalciumOnly);
    abn_idx    = find(isAbnormalFit);
    plot_order = [norm_idx(:)', gray_idx(:)', ca_idx(:)', abn_idx(:)'];

    % Colori assi per categoria
    col_norm = [0   0   0  ];   % nero
    col_gray = [0.5 0.5 0.5];   % grigio
    col_ca   = [0   0   0.8];   % blu
    col_abn  = [0.8 0   0  ];   % rosso

    figTitle = sprintf('Riga %i | Mouse: %s | Field: %s | Stim: %s | Scaling: %s', ...
        iRow, ClCa(iRow).Mouse, ClCa(iRow).Field, ClCa(iRow).Stimulus, scalingMethod);

    figure('Name', figTitle, 'NumberTitle', 'off', 'Color', 'w');
    tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

    for roi = plot_order
        nexttile

        plot(Time_920, mean920_all(:, roi),   'k', 'LineWidth', 1.2); hold on
        plot(Time_920, scaledFit_all(:, roi), 'r', 'LineWidth', 1.5);
        plot(Time_920, residual_all(:, roi),  'b', 'LineWidth', 1.2);
        xline(0, 'k--'); yline(0, 'k:');
        xlim([min(Time_920) max(Time_920)])
        set(gca, 'Box', 'off')

        if isNormal(roi)
            col = col_norm; tag = '';
        elseif isGray(roi)
            col = col_gray; tag = ' (Gray)';
        elseif isCalciumOnly(roi)
            col = col_ca;   tag = ' (Ca-only)';
        else
            col = col_abn;  tag = ' (ABNORMAL)';
        end

        set(gca, 'XColor', col, 'YColor', col)
        title(sprintf('ROI %i%s  \\alpha=%.2f', roi, tag, scaleFactor_all(roi)), ...
            'Color', col, 'FontSize', 7)

        if roi == plot_order(1)
            legend({'920nm dG/R','Scaled avatar 830','Residual'}, ...
                'Location','best','FontSize',7);
        end
    end
    sgtitle(figTitle, 'Interpreter', 'none')

end % fine loop righe

fprintf('--- Elaborazione completata ---\n\n');


% =========================================================================
%% SEZIONE 3 — SALVATAGGIO
% =========================================================================

risposta_salva = questdlg('Vuoi salvare ClCa aggiornata?', ...
    'Salvataggio', 'Sì', 'No', 'Sì');
if strcmp(risposta_salva, 'Sì')
    defaultName = ['ClCa_session3_', strrep(strjoin(stimToAnalyze,'+'),'_',''), '_', ...
                   datestr(now, 'yyyymmdd_HHMM'), '.mat']; %#ok<TNOW1,DATST>
    [saveName, saveDir] = uiputfile('*.mat', 'Salva ClCa aggiornata', defaultName);
    if ~isequal(saveName, 0)
        save(fullfile(saveDir, saveName), 'ClCa', '-v7.3');
        fprintf('ClCa salvata in:\n  %s\n', fullfile(saveDir, saveName));
    else
        fprintf('Salvataggio annullato.\n');
    end
else
    fprintf('ClCa disponibile nel workspace.\n');
end

fprintf('\n=========================================================\n');
fprintf('               Sessione 3 completata.                    \n');
fprintf('=========================================================\n');
fprintf('\nVariabili disponibili per la sessione 4:\n');
fprintf('  stimToAnalyze = {%s}\n', strjoin(stimToAnalyze, ', '));
fprintf('  rowsWithStim  = [%s]\n', num2str(rowsWithStim));


%% =========================================================================
% =========================================================================
% ClCa_session4_GCaMPcorrection.m  (G/R version)
% =========================================================================

fprintf('=========================================================\n');
fprintf('   ClCa SESSION 4 - Correzione GCaMP e deltaG/G  (G/R)   \n');
fprintf('=========================================================\n\n');

if ~exist('ClCa','var') || isempty(ClCa)
    error('ClCa non trovata. Esegui prima la sessione 3.');
end
if ~exist('rowsWithStim','var') || isempty(rowsWithStim)
    error('rowsWithStim non trovata. Esegui prima la sessione 3.');
end
if ~exist('stimToAnalyze','var') || isempty(stimToAnalyze)
    error('stimToAnalyze non trovata. Esegui prima la sessione 3.');
end

idxRows = rowsWithStim;


% =========================================================================
%% SEZIONE 0 — PARAMETRI EDITABILI
% =========================================================================

K830 = 1.4;             % <-- MODIFICA QUI
K920 = 0.85;            % <-- MODIFICA QUI
alpha_bleed = 0;        % <-- MODIFICA QUI  bleedthrough (rinominato per
%                                            evitare conflitto con alpha=scaleFactor)
mov_avg_win = [10, 0];  % <-- MODIFICA QUI  smooth figure intermedie

% Finestra post-stimolo per SNR [inizio, fine] in secondi
snr_signal_window = [0, 4.5];   % <-- MODIFICA QUI

fprintf('K830=%.3f | K920=%.3f | bleedthrough=%.3f\n', K830, K920, alpha_bleed);
fprintf('Finestra SNR: [%.1f, %.1f] s\n\n', snr_signal_window(1), snr_signal_window(2));

% Colori categoria (coerenti con session3)
col_norm = [0   0   0  ];
col_gray = [0.5 0.5 0.5];
col_ca   = [0   0   0.8];
col_abn  = [0.8 0   0  ];


% =========================================================================
% =========================================================================
%% SEZIONE 1 — CORREZIONE E CALCOLO deltaG/G
% =========================================================================

fprintf('--- Sezione 1: Correzione GreenGCaMP e deltaG/G ---\n\n');

for k = 1 : numel(idxRows)

    i           = idxRows(k);
    fieldSuffix = strrep(ClCa(i).Stimulus, '_', '');
    res         = ['Results_' fieldSuffix];

    % Verifica che il campo Results esista per questa riga
    if ~isfield(ClCa(i), res)
        fprintf('  [SKIP] Riga %i: campo %s non trovato. Riesegui session3.\n\n', i, res);
        continue
    end

    % Vettori tempo separati per i due canali
    % (possono avere frequenze diverse → nTime diverso)
    Time_830 = ClCa(i).(res).Time_830(:);
    Time_920 = ClCa(i).(res).Time_920(:);

    baselineIdx_830 = Time_830 < 0;
    baselineIdx_920 = Time_920 < 0;

    G830 = ClCa(i).G_smooth_830;   % nTime_830 x nTrials_830 x nROIs
    R830 = ClCa(i).R_smooth_830;
    G920 = ClCa(i).G_smooth_920;   % nTime_920 x nTrials_920 x nROIs
    R920 = ClCa(i).R_smooth_920;

    nROIs       = size(G830, 3);
    nTrials_830 = size(G830, 2);
    nTrials_920 = size(G920, 2);

    Avatar        = ClCa(i).(res).ScaledFit;   % nTime_920 x nROIs
    isNormal      = ClCa(i).(res).IsNormal;
    isCalciumOnly = ClCa(i).(res).IsCalciumOnly;
    isAbnormalFit = ClCa(i).(res).IsAbnormalFit;
    isGray        = ClCa(i).(res).IsGray;
    noCorrection  = ~isNormal;

    % ---- Media sui trial per ciascun canale -----------------------------
    % Matematicamente equivalente a correggere trial per trial e poi mediare
    % perché la correzione è lineare. Robusto a nTrials diversi tra canali.

    G920_mean = squeeze(mean(G920, 2, 'omitnan'));  % nTime_920 x nROIs
    R920_mean = squeeze(mean(R920, 2, 'omitnan'));
    G830_mean = squeeze(mean(G830, 2, 'omitnan'));  % nTime_830 x nROIs
    R830_mean = squeeze(mean(R830, 2, 'omitnan'));

    % Gestisce il caso nROIs=1 (squeeze riduce a vettore colonna)
    if nROIs == 1
        G920_mean = G920_mean(:);  R920_mean = R920_mean(:);
        G830_mean = G830_mean(:);  R830_mean = R830_mean(:);
    end

    % ---- Correzione sulla media ------------------------------------------
    % Avatar (nTime_920 x nROIs) è compatibile con G920_mean/R920_mean.
    % Per 830nm, l'Avatar va interpolato su Time_830 perché ScaledFit
    % è costruito sul tempo di 920nm.
    Avatar_on830 = interp1(Time_920, Avatar, Time_830, 'linear', 'extrap');

    GCaMP_920_mean = G920_mean - (K920 .* (1 + Avatar)      + alpha_bleed) .* R920_mean;
    GCaMP_830_mean = G830_mean - (K830 .* (1 + Avatar_on830) + alpha_bleed) .* R830_mean;

    % ROI non-normali: nessuna correzione (GCaMP = G grezzo mediato)
    for r = 1:nROIs
        if noCorrection(r)
            GCaMP_920_mean(:, r) = G920_mean(:, r);
            GCaMP_830_mean(:, r) = G830_mean(:, r);
        end
    end

    % ---- Salva versione per trial (utile per SEM nelle figure) ----------
    Avatar_920 = repmat(permute(Avatar,      [1 3 2]), [1 nTrials_920 1]);
    Avatar_830 = repmat(permute(Avatar_on830,[1 3 2]), [1 nTrials_830 1]);

    GCaMP_920_trials = G920 - (K920 .* (1 + Avatar_920) + alpha_bleed) .* R920;
    GCaMP_830_trials = G830 - (K830 .* (1 + Avatar_830) + alpha_bleed) .* R830;
    for r = 1:nROIs
        if noCorrection(r)
            GCaMP_920_trials(:,:,r) = G920(:,:,r);
            GCaMP_830_trials(:,:,r) = G830(:,:,r);
        end
    end

    ClCa(i).GreenGCaMP_920 = GCaMP_920_trials;
    ClCa(i).GreenGCaMP_830 = GCaMP_830_trials;

    % ---- Calcolo deltaG/G -----------------------------------------------
    G0_920 = mean(GCaMP_920_mean(baselineIdx_920, :), 1, 'omitnan');  % 1 x nROIs
    G0_830 = mean(GCaMP_830_mean(baselineIdx_830, :), 1, 'omitnan');

    ClCa(i).DeltaGG_920 = 100 .* (GCaMP_920_mean - G0_920) ./ G0_920;
    ClCa(i).DeltaGG_830 = 100 .* (GCaMP_830_mean - G0_830) ./ G0_830;

    fprintf('  Riga %i | Normali:%i | Ca-only:%i | Abnormal:%i | Grigie:%i\n', ...
        i, sum(isNormal), sum(isCalciumOnly), sum(isAbnormalFit), sum(isGray));

end
fprintf('\n');


% =========================================================================
%% SEZIONE 2 — FIGURE INTERMEDIE (opzionali)
% =========================================================================

if strcmp(plotIntermediate, 'Sì')
    for k = 1 : numel(idxRows)
        i           = idxRows(k);
        fieldSuffix = strrep(ClCa(i).Stimulus, '_', '');
        res         = ['Results_' fieldSuffix];
        if ~isfield(ClCa(i), res), continue; end

        Time_830 = ClCa(i).(res).Time_830(:);
        Time_920 = ClCa(i).(res).Time_920(:);
        nROIs       = size(ClCa(i).GreenGCaMP_920, 3);
        nTrials_920 = size(ClCa(i).GreenGCaMP_920, 2);
        nTrials_830 = size(ClCa(i).GreenGCaMP_830, 2);
        isNormal      = ClCa(i).(res).IsNormal;
        isCalciumOnly = ClCa(i).(res).IsCalciumOnly;
        isAbnormalFit = ClCa(i).(res).IsAbnormalFit;
        isGray        = ClCa(i).(res).IsGray;

        norm_idx   = find(isNormal);
        gray_idx   = find(isGray);
        ca_idx     = find(isCalciumOnly);
        abn_idx    = find(isAbnormalFit);
        plot_order = [norm_idx(:)', gray_idx(:)', ca_idx(:)', abn_idx(:)'];

        rowID = sprintf('Riga %i | Mouse:%s | Field:%s | Stim:%s', ...
            i, ClCa(i).Mouse, ClCa(i).Field, ClCa(i).Stimulus);

        % G e R smoothati 920nm (usa Time_920)
        figure('Name', ['G/R smooth 920nm | ' rowID], 'Color','w');
        tiledlayout('flow','TileSpacing','compact');
        for r = plot_order
            nexttile; hold on
            g_m = movmean(mean(ClCa(i).G_smooth_920(:,:,r), 2,'omitnan'), mov_avg_win);
            r_m = movmean(mean(ClCa(i).R_smooth_920(:,:,r), 2,'omitnan'), mov_avg_win);
            plot(Time_920, g_m, 'g', 'LineWidth', 1.5)
            plot(Time_920, r_m, 'r', 'LineWidth', 1.5)
            xline(0,'k--'); axis tight; set(gca,'Box','off');
            if isNormal(r),          col = col_norm; tag = '';
            elseif isGray(r),        col = col_gray; tag = ' (Gray)';
            elseif isCalciumOnly(r), col = col_ca;   tag = ' (Ca)';
            else,                    col = col_abn;  tag = ' (Abn)';
            end
            set(gca, 'XColor', col, 'YColor', col);
            title(sprintf('ROI %i%s', r, tag), 'Color', col, 'FontSize', 7);
        end
        sgtitle(['920nm G & R | ' rowID], 'Interpreter','none');
    end
end


% =========================================================================
% =========================================================================
%% SEZIONE 3 — SNR E FIGURA ORDINATA (sempre eseguita)
% =========================================================================

fprintf('--- Calcolo SNR e figura ordinata ---\n\n');

for k = 1 : numel(idxRows)

    i           = idxRows(k);
    fieldSuffix = strrep(ClCa(i).Stimulus, '_', '');
    res         = ['Results_' fieldSuffix];

    if ~isfield(ClCa(i), res), continue; end

    % Vettori tempo separati — regola generale:
    % DeltaGG_920 va sempre accoppiato con Time_920
    % DeltaGG_830 va sempre accoppiato con Time_830
    Time_830 = ClCa(i).(res).Time_830(:);
    Time_920 = ClCa(i).(res).Time_920(:);

    isNormal      = ClCa(i).(res).IsNormal;
    isCalciumOnly = ClCa(i).(res).IsCalciumOnly;
    isAbnormalFit = ClCa(i).(res).IsAbnormalFit;
    isGray        = ClCa(i).(res).IsGray;
    nROIs         = size(ClCa(i).DeltaGG_920, 2);

    % Maschere temporali su Time_920 (per DeltaGG_920 usato nell'SNR)
    baseline_mask_920 = (Time_920 < 0) | (Time_920 > 20);
    signal_mask_920   = (Time_920 >= snr_signal_window(1)) & ...
                        (Time_920 <= snr_signal_window(2));

    % Calcolo SNR per ogni ROI su DeltaGG_920
    snr_all = nan(nROIs, 1);
    for r = 1:nROIs
        dgg        = ClCa(i).DeltaGG_920(:, r);   % nTime_920 x 1
        noise_val  = std(dgg(baseline_mask_920), 'omitnan');
        signal_val = abs(mean(dgg(signal_mask_920), 'omitnan'));
        if noise_val > 0 && ~isnan(noise_val)
            snr_all(r) = signal_val / noise_val;
        end
    end
    ClCa(i).(res).SNR_920 = snr_all;

    % Ordina ogni gruppo per SNR decrescente
    norm_idx = find(isNormal);
    gray_idx = find(isGray);
    ca_idx   = find(isCalciumOnly);
    abn_idx  = find(isAbnormalFit);

    sort_group = @(idx) idx(argsort_desc(snr_all(idx)));

    sorted_norm = sort_group(norm_idx);
    sorted_gray = sort_group(gray_idx);
    sorted_ca   = sort_group(ca_idx);
    sorted_abn  = sort_group(abn_idx);

    sorted_order = [sorted_norm(:)', sorted_gray(:)', ...
                    sorted_ca(:)', sorted_abn(:)'];

    % ---- Figura ordinata ------------------------------------------------
    rowID = sprintf('Riga %i | Mouse:%s | Field:%s | Stim:%s', ...
        i, ClCa(i).Mouse, ClCa(i).Field, ClCa(i).Stimulus);

    figure('Name', ['DeltaGG SORTED by SNR | ' rowID], 'Color','w');
    tiledlayout('flow','TileSpacing','compact');

    for step = 1:length(sorted_order)
        r = sorted_order(step);
        nexttile; hold on

        % DeltaGG_920 → Time_920 | DeltaGG_830 → Time_830
        plot(Time_920, ClCa(i).DeltaGG_920(:,r), ...
            'Color',[0.0 0.6 0.0], 'LineWidth',1.5);
        plot(Time_830, ClCa(i).DeltaGG_830(:,r), ...
            'Color',[0.5 0.9 0.5], 'LineWidth',1.2);
        xline(0,'k--'); yline(0,'k:');

        % Usa il range temporale comune per xlim
        t_start = min(Time_920(1),   Time_830(1));
        t_end   = max(Time_920(end), Time_830(end));
        xlim([t_start, t_end]);
        set(gca,'Box','off');

        if isNormal(r),          col = col_norm; tag = '';
        elseif isGray(r),        col = col_gray; tag = ' Gray';
        elseif isCalciumOnly(r), col = col_ca;   tag = ' Ca-only';
        else,                    col = col_abn;  tag = ' Abnormal';
        end

        title(sprintf('ROI %i%s (SNR:%.2f)', r, tag, snr_all(r)), ...
            'Color', col, 'FontSize', 8);
        set(gca, 'XColor', col, 'YColor', col);

        if step == 1
            legend({'920nm','830nm'}, 'FontSize',6, 'Location','best');
        end
        ylabel('\DeltaG/G (%)');
    end
    sgtitle(['\DeltaG/G per SNR | ' rowID], 'Interpreter','none');

    fprintf('  Riga %i | Normali:%i | Grigie:%i | Ca-only:%i | Abnormal:%i\n', ...
        i, numel(sorted_norm), numel(sorted_gray), ...
        numel(sorted_ca), numel(sorted_abn));
end
fprintf('\n');

% =========================================================================
%% SEZIONE 4 — SALVATAGGIO
% =========================================================================

risposta_salva = questdlg('Vuoi salvare ClCa aggiornata?', 'Salvataggio', 'Sì', 'No', 'Sì');
if strcmp(risposta_salva, 'Sì')
    defaultName = sprintf('ClCa_%s_session4_%s.mat', ...
        strjoin(stimToAnalyze,'+'), datestr(now,'yyyymmdd_HHMM')); %#ok<TNOW1,DATST>
    [saveName, saveDir] = uiputfile('*.mat','Salva ClCa', defaultName);
    if ~isequal(saveName, 0)
        save(fullfile(saveDir, saveName), 'ClCa', '-v7.3');
        fprintf('ClCa salvata in:\n  %s\n', fullfile(saveDir, saveName));
    else
        fprintf('Salvataggio annullato.\n');
    end
else
    fprintf('ClCa disponibile nel workspace.\n');
end

fprintf('\n=========================================================\n');
fprintf('               Sessione 4 completata.                    \n');
fprintf('=========================================================\n');


% =========================================================================





%% =========================================================================
% ClCa_session5_unmaskedAnalysis.m
% =========================================================================
%
% DESCRIZIONE:
%   Quinta sessione di analisi sulla struttura ClCa.
%   Va eseguita dopo la sessione 4 (correzione GCaMP e deltaG/G).
%
%   OPERAZIONI:
%
%   SEZIONE 1 — IDENTIFICAZIONE ROI UNMASKED
%     Per ogni riga con lo stimolo analizzato, tra le ROI "normali"
%     (asse nero), identifica quelle che:
%       a) SNR > snr_threshold         (risposta affidabile)
%       b) Picco DeltaGG_920 > calcium_response_threshold
%          nella finestra [0, calcium_search_end_sec]   (risposta calcio)
%     Queste ROI sono dette "Unmasked": la correzione con l'avatar ha
%     smascherato la loro risposta calcio che era nascosta dal segnale cloro.
%
%   SEZIONE 2 — PLOT ROI UNMASKED (opzionale)
%     Plotta tutte le ROI unmasked di tutte le righe, 70 per figura.
%     Titolo di ogni tile: riga e numero ROI.
%
%   SEZIONE 3 — DEBUG SINGOLA ROI (opzionale)
%     Per una riga e ROI specificate dall'operatore, mostra:
%       - Fig A: dG/R 830nm, dG/R 920nm, avatar scalato
%       - Fig B: DeltaGG 830nm e 920nm
%
%   SEZIONE 4 — FIGURE RIASSUNTIVE
%     Fig 1: media ± SD delle ROI unmasked
%       - Subplot 1: DeltaGG_920 corretto (GCaMP post-correzione)
%       - Subplot 2: DeltaGG_920 raw (GCaMP senza correzione avatar,
%                    ma già calcolato in session4 per le Ca-only)
%     Fig 2: media ± SD delle ROI solo-cloro (IsNormal ma non unmasked,
%             con risposta cloro ma senza risposta calcio smascherata)
%             usando DeltaGG_920 raw (nessuna correzione applicata a queste)
%
% PARAMETRI EDITABILI (Sezione 0):
%   snr_threshold            : SNR minimo per considerare una ROI unmasked
%   calcium_search_end_sec   : fine finestra ricerca risposta calcio [s]
%   calcium_response_threshold: picco minimo DeltaGG_920 per risposta calcio [%]
%
% INPUT:
%   ClCa con i campi prodotti da session3 e session4.
%   rowsWithStim e stimToAnalyze dal workspace (session3).
%
% =========================================================================

fprintf('=========================================================\n');
fprintf('   ClCa SESSION 5 - ROI Unmasked e Figure Riassuntive    \n');
fprintf('=========================================================\n\n');

% Verifica variabili workspace
if ~exist('ClCa','var') || isempty(ClCa)
    error('ClCa non trovata. Esegui prima le sessioni 3 e 4.');
end
if ~exist('rowsWithStim','var') || isempty(rowsWithStim)
    error('rowsWithStim non trovata. Esegui prima la sessione 3.');
end
if ~exist('stimToAnalyze','var') || isempty(stimToAnalyze)
    error('stimToAnalyze non trovata. Esegui prima la sessione 3.');
end

idxRows = rowsWithStim;
fprintf('Stimoli: %s\n', strjoin(stimToAnalyze, ', '));
fprintf('Righe da analizzare: %i\n\n', numel(idxRows));

% Colori categoria (coerenti con session3/4)
col_norm = [0   0   0  ];
col_gray = [0.5 0.5 0.5];
col_ca   = [0   0   0.8];
col_abn  = [0.8 0   0  ];
col_unm  = [0.1 0.7 0.1];  % verde scuro per unmasked


% =========================================================================
%% SEZIONE 0 — PARAMETRI EDITABILI
% =========================================================================

% SNR minimo (calcolato su DeltaGG_920) per considerare una ROI unmasked.
% Solo le ROI normali con SNR > soglia vengono candidate.
snr_threshold = 1.5;             % <-- MODIFICA QUI

% Fine della finestra temporale in cui cercare il picco calcio [secondi].
% La risposta GCaMP è veloce: tipicamente tra 0 e 3-4 secondi.
calcium_search_end_sec = 3.5;    % <-- MODIFICA QUI

% Ampiezza minima del picco DeltaGG_920 nella finestra post-stimolo
% per considerare la risposta calcio significativa [%].
% Es. 10 = picco di almeno il 10% sopra la baseline.
calcium_response_threshold = 20; % <-- MODIFICA QUI  [%]

% Soglia massima ammissibile per DeltaG/G. 
% Se 1.0 = 100%, un valore di 5.0 significa 500%. Oltre è quasi certamente un artefatto.
max_deltaGG_limit = 2000.0; 

% Soglia massima per la Deviazione Standard.
% Una ROI "sana" ha un rumore di fondo limitato. Se la STD è enorme, la traccia è impazzita.
max_std_limit = 20000.0;

fprintf('Parametri identificazione Unmasked:\n');
fprintf('  SNR minimo                  : %.2f\n', snr_threshold);
fprintf('  Finestra ricerca calcio     : [0, %.1f] s\n', calcium_search_end_sec);
fprintf('  Soglia risposta calcio      : %.1f %%\n\n', calcium_response_threshold);
fprintf('  Soglia sanità      : %.2f\n', max_deltaGG_limit);


% =========================================================================
%% SEZIONE 1 — IDENTIFICAZIONE ROI UNMASKED
% =========================================================================

fprintf('--- Sezione 1: Identificazione ROI Unmasked ---\n\n');

for k = 1 : numel(idxRows)

    i           = idxRows(k);
    fieldSuffix = strrep(ClCa(i).Stimulus, '_', '');
    res         = ['Results_' fieldSuffix];
    if ~isfield(ClCa(i), res), continue; end

    % Time_920 perché IsUnmasked si basa su DeltaGG_920
    Time_920   = ClCa(i).(res).Time_920(:);
    isNormal   = ClCa(i).(res).IsNormal;
    snr_all    = ClCa(i).(res).SNR_920;
    nROIs      = size(ClCa(i).DeltaGG_920, 2);

    calciumWin = (Time_920 >= 0) & (Time_920 <= calcium_search_end_sec);

    isUnmasked = false(nROIs, 1);

    for r = 1:nROIs
        if ~isNormal(r), continue; end
        if isnan(snr_all(r)) || snr_all(r) <= snr_threshold, continue; end

        dgg920_roi   = ClCa(i).DeltaGG_920(:, r);   % nTime_920 x 1
        % Calcoliamo i parametri di "salute" della traccia
        current_max = max(abs(dgg920_roi), [], 'omitnan');
        current_std = std(dgg920_roi, 0, 'omitnan');

        % Verifichiamo se la ROI è "fuori controllo"
        is_anomalous = (current_max > max_deltaGG_limit) || (current_std > max_std_limit);

        if is_anomalous
            % Se la ROI è anomala, la scartiamo a prescindere e passiamo alla prossima
            % (Opzionale: puoi loggare quale ROI è stata scartata per anomalie)
            % fprintf('Nota: ROI %d riga %d scartata per valori fuori scala (Max:%.1f, STD:%.1f)\n', r, i, current_max, current_std);
            continue;
        end
        % --------------------------------------------

        % Se passa il filtro sopra, procediamo con il controllo del calcio
        peak_calcium = max(dgg920_roi(calciumWin));

        if peak_calcium >= calcium_response_threshold
            isUnmasked(r) = true;
        end
    end

    ClCa(i).(res).IsUnmasked = isUnmasked;

    fprintf('  Riga %i | Mouse:%s | Field:%s | Stim:%s\n', ...
        i, ClCa(i).Mouse, ClCa(i).Field, ClCa(i).Stimulus);
    fprintf('    ROI normali: %i | ROI unmasked: %i\n\n', ...
        sum(isNormal), sum(isUnmasked));
end

fprintf('Identificazione completata.\n\n');


% =========================================================================
%% SEZIONE 2 — PLOT ROI UNMASKED (opzionale, 70 per figura)
% =========================================================================

plotUnmasked = questdlg('Vuoi plottare tutte le ROI unmasked?', ...
    'Plot ROI Unmasked', 'Sì', 'No', 'No');

if strcmp(plotUnmasked, 'Sì')

    fprintf('--- Sezione 2: Plot ROI Unmasked ---\n');

    % Raccoglie tutte le ROI unmasked in ordine: per riga, per roi
    unmasked_list = [];  % struttura: [iRow, roi]
    for k = 1 : numel(idxRows)
        i           = idxRows(k);
        fieldSuffix = strrep(ClCa(i).Stimulus, '_', '');
        res         = ['Results_' fieldSuffix];
        isUnmasked  = ClCa(i).(res).IsUnmasked;
        roi_list    = find(isUnmasked);
        for r = roi_list(:)'
            unmasked_list = [unmasked_list; i, r]; %#ok<AGROW>
        end
    end

    nUnmasked   = size(unmasked_list, 1);
    rois_per_fig = 70;
    nFigs        = ceil(nUnmasked / rois_per_fig);

    fprintf('  ROI unmasked totali: %i | Figure necessarie: %i\n\n', nUnmasked, nFigs);

    for fig_num = 1 : nFigs
        idx_start = (fig_num - 1) * rois_per_fig + 1;
        idx_end   = min(fig_num * rois_per_fig, nUnmasked);
        batch     = unmasked_list(idx_start:idx_end, :);

        figure('Name', sprintf('ROI Unmasked | Fig %i/%i', fig_num, nFigs), ...
               'NumberTitle','off','Color','w');
        tiledlayout('flow','TileSpacing','compact','Padding','compact');

        for j = 1 : size(batch, 1)
    iRow = batch(j, 1);
    r    = batch(j, 2);

    fieldSuffix = strrep(ClCa(iRow).Stimulus, '_', '');
    res         = ['Results_' fieldSuffix];
    % DeltaGG_920 → Time_920 | DeltaGG_830 → Time_830
    Time_920_p  = ClCa(iRow).(res).Time_920(:);
    Time_830_p  = ClCa(iRow).(res).Time_830(:);
    snr_val     = ClCa(iRow).(res).SNR_920(r);

    nexttile; hold on
    plot(Time_920_p, ClCa(iRow).DeltaGG_920(:, r), ...
        'Color', col_unm, 'LineWidth', 1.5);
    plot(Time_830_p, ClCa(iRow).DeltaGG_830(:, r), ...
        'Color', [0.5 0.9 0.5], 'LineWidth', 1.0);
    xline(0,'k--'); yline(0,'k:');
    t_start = min(Time_920_p(1),   Time_830_p(1));
    t_end   = max(Time_920_p(end), Time_830_p(end));
    xlim([t_start, t_end]);
    set(gca,'Box','off');

    title(sprintf('R%i ROI%i (SNR:%.1f)', iRow, r, snr_val), ...
        'FontSize', 7, 'Color', col_unm);
    ylabel('\DeltaG/G (%)');
end

        sgtitle(sprintf('ROI Unmasked | Stim: %s | Fig %i/%i', ...
            strjoin(stimToAnalyze,'+'), fig_num, nFigs), ...
            'Interpreter','none');
    end

    fprintf('  Plot completato.\n\n');
end


% =========================================================================
%% SEZIONE 3 — DEBUG SINGOLA ROI (opzionale)
% =========================================================================

debugROI = questdlg('Vuoi eseguire il debug di una singola ROI?', ...
    'Debug ROI', 'Sì', 'No', 'No');

if strcmp(debugROI, 'Sì')

    fprintf('--- Sezione 3: Debug singola ROI ---\n');

    % Chiede riga e ROI all'operatore
    answer = inputdlg({'Numero riga (indice ClCa):', 'Numero ROI:'}, ...
        'Debug ROI', 1, {'1','1'});
    if isempty(answer)
        fprintf('  Debug annullato.\n\n');
    else
        iRow_dbg = str2double(answer{1});
        roi_dbg  = str2double(answer{2});

        if isnan(iRow_dbg) || isnan(roi_dbg)
            fprintf('  [ERRORE] Valori non validi inseriti.\n\n');
        else
            fieldSuffix = strrep(ClCa(iRow_dbg).Stimulus, '_', '');
            res         = ['Results_' fieldSuffix];

            Time_830 = ClCa(iRow_dbg).(res).Time_830(:);
            Time_920 = ClCa(iRow_dbg).(res).Time_920(:);

            dRatio_830 = ClCa(iRow_dbg).(res).dRatio_830(:, roi_dbg);
            dRatio_920 = ClCa(iRow_dbg).(res).dRatio_920(:, roi_dbg);
            scaledFit  = ClCa(iRow_dbg).(res).ScaledFit(:, roi_dbg);
            sf_val     = ClCa(iRow_dbg).(res).ScaleFactor(roi_dbg);
            snr_val    = ClCa(iRow_dbg).(res).SNR_920(roi_dbg);

            % Determina categoria
            if ClCa(iRow_dbg).(res).IsUnmasked(roi_dbg)
                cat_str = 'UNMASKED'; cat_col = col_unm;
            elseif ClCa(iRow_dbg).(res).IsNormal(roi_dbg)
                cat_str = 'Normale'; cat_col = col_norm;
            elseif ClCa(iRow_dbg).(res).IsCalciumOnly(roi_dbg)
                cat_str = 'Ca-only'; cat_col = col_ca;
            elseif ClCa(iRow_dbg).(res).IsAbnormalFit(roi_dbg)
                cat_str = 'Abnormal'; cat_col = col_abn;
            else
                cat_str = 'Gray'; cat_col = col_gray;
            end

            dbgTitle = sprintf('Riga %i | ROI %i | %s | \\alpha=%.2f | SNR=%.2f', ...
                iRow_dbg, roi_dbg, cat_str, sf_val, snr_val);

            % --- Figura A: dG/R e avatar scalato -------------------------
            figure('Name', ['DEBUG Fig A | ' dbgTitle], ...
                   'NumberTitle','off','Color','w');
            tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');

            nexttile; hold on
            plot(Time_830, dRatio_830, 'k', 'LineWidth', 1.5);
            plot(Time_920, scaledFit,  'r', 'LineWidth', 2.0);
            xline(0,'k--'); yline(0,'k:');
            legend({'dG/R 830nm (raw mean)','Avatar scalato 830→920'}, ...
                'Location','best','FontSize',8);
            xlabel('Time (s)'); ylabel('dG/R (normalizzato)');
            title('830nm: segnale iClima + avatar');
            set(gca,'Box','off'); xlim([Time_830(1), Time_830(end)]);

            nexttile; hold on
            plot(Time_920, dRatio_920, 'Color',[0 0 0.6], 'LineWidth', 1.5);
            plot(Time_920, scaledFit,  'r', 'LineWidth', 2.0);
            xline(0,'k--'); yline(0,'k:');
            legend({'dG/R 920nm (raw mean)','Avatar scalato 830→920'}, ...
                'Location','best','FontSize',8);
            xlabel('Time (s)'); ylabel('dG/R (normalizzato)');
            title('920nm: segnale misto + avatar');
            set(gca,'Box','off'); xlim([Time_920(1), Time_920(end)]);

            sgtitle(['Fig A - dG/R e Avatar | ' dbgTitle], ...
                'Interpreter','none','Color',cat_col);

            % --- Figura B: DeltaGG 830 e 920 -----------------------------
            figure('Name', ['DEBUG Fig B | ' dbgTitle], ...
       'NumberTitle','off','Color','w');
hold on
% DeltaGG_920 → Time_920 | DeltaGG_830 → Time_830
Time_920_d = ClCa(iRow_dbg).(res).Time_920(:);
Time_830_d = ClCa(iRow_dbg).(res).Time_830(:);

plot(Time_920_d, ClCa(iRow_dbg).DeltaGG_920(:, roi_dbg), ...
    'Color',[0.0 0.6 0.0], 'LineWidth', 2.0);
plot(Time_830_d, ClCa(iRow_dbg).DeltaGG_830(:, roi_dbg), ...
    'Color',[0.5 0.9 0.5], 'LineWidth', 1.5);
xline(0,'k--'); yline(0,'k:');
legend({'\DeltaG/G 920nm (corretto)','\DeltaG/G 830nm (corretto)'}, ...
    'Location','best','FontSize',9);
xlabel('Time (s)'); ylabel('\DeltaG/G (%)');
title(['Fig B - \DeltaG/G | ' dbgTitle], 'Interpreter','none','Color',cat_col);
set(gca,'Box','off');
t_start = min(Time_920_d(1),   Time_830_d(1));
t_end   = max(Time_920_d(end), Time_830_d(end));
xlim([t_start, t_end]);

            fprintf('  Debug ROI %i di riga %i completato.\n\n', roi_dbg, iRow_dbg);
        end
    end
end


% =========================================================================
%% SEZIONE 4 — FIGURE RIASSUNTIVE
% =========================================================================

fprintf('--- Sezione 4: Figure riassuntive ---\n\n');

% ---- Raccoglie tutte le ROI unmasked (pool globale) ---------------------
% Ogni ROI contribuisce con la sua traccia DeltaGG_920 corretta e raw.
% Le tracce vengono interpolate su un asse temporale comune (Time_ref)
% per permettere la media tra righe con frequenze diverse.

% Costruisce Time_ref come unione degli assi temporali di tutte le righe
% (in pratica: prende il range comune e usa la risoluzione più fine)
t_starts = zeros(numel(idxRows), 1);
t_ends   = zeros(numel(idxRows), 1);
dt_all   = zeros(numel(idxRows), 1);

for k = 1:numel(idxRows)
    i           = idxRows(k);
    fieldSuffix = strrep(ClCa(i).Stimulus, '_', '');
    res         = ['Results_' fieldSuffix];
    if ~isfield(ClCa(i), res), continue; end
    % Usa Time_920 per il pool perché DeltaGG_920 è il segnale principale
    T            = ClCa(i).(res).Time_920(:);
    t_starts(k)  = T(1);
    t_ends(k)    = T(end);
    dt_all(k)    = median(diff(T));
end

t_min    = max(t_starts);
t_max    = min(t_ends);
dt_ref   = min(dt_all(dt_all > 0));
Time_ref = (t_min : dt_ref : t_max)';

pool_unmasked_corr = [];
pool_unmasked_raw  = [];
pool_calcio_raw    = [];
pool_unmasked_GR830 = [];
pool_unmasked_GR920 = [];

for k = 1 : numel(idxRows)
    i           = idxRows(k);
    fieldSuffix = strrep(ClCa(i).Stimulus, '_', '');
    res         = ['Results_' fieldSuffix];
    if ~isfield(ClCa(i), res), continue; end

    Time_920_i    = ClCa(i).(res).Time_920(:);
    Time_830_i    = ClCa(i).(res).Time_830(:);
    isUnmasked    = ClCa(i).(res).IsUnmasked;
    isNormal      = ClCa(i).(res).IsNormal;
    isCalciumOnly = ClCa(i).(res).IsCalciumOnly;

    % Pool DeltaGG_920 corretto (usa Time_920)
    for r = find(isUnmasked)'
        tr = interp1(Time_920_i, ClCa(i).DeltaGG_920(:,r), Time_ref, 'linear', NaN);
        pool_unmasked_corr = [pool_unmasked_corr, tr]; %#ok<AGROW>
    end

    % Pool raw (G_smooth_920 non corretto, normalizzato)
    for r = find(isUnmasked)'
        G920_r    = squeeze(ClCa(i).G_smooth_920(:,:,r));
        G920_mean = mean(G920_r, 2, 'omitnan');
        baseIdx   = Time_920_i < 0;
        G0        = mean(G920_mean(baseIdx), 'omitnan');
        dgg_raw   = 100 .* (G920_mean - G0) ./ G0;
        tr        = interp1(Time_920_i, dgg_raw, Time_ref, 'linear', NaN);
        pool_unmasked_raw = [pool_unmasked_raw, tr]; %#ok<AGROW>
    end

    % Pool solo-calcio (IsCalciumOnly, DeltaGG_920, usa Time_920)
    for r = find(isCalciumOnly)'
        tr = interp1(Time_920_i, ClCa(i).DeltaGG_920(:,r), Time_ref, 'linear', NaN);
        pool_calcio_raw = [pool_calcio_raw, tr]; %#ok<AGROW>
    end

    for r = find(isUnmasked)'
        % 830 nm (usa dRatio_830 e il suo asse temporale specifico)
        tr830 = interp1(Time_830_i, ClCa(i).(res).dRatio_830(:,r), Time_ref, 'linear', NaN);
        pool_unmasked_GR830 = [pool_unmasked_GR830, tr830 .* 100]; % Moltiplico per 100 per avere %
        
        % 920 nm (usa dRatio_920 e il suo asse temporale specifico)
        tr920 = interp1(Time_920_i, ClCa(i).(res).dRatio_920(:,r), Time_ref, 'linear', NaN);
        pool_unmasked_GR920 = [pool_unmasked_GR920, tr920 .* 100];
    end
end

nUnmasked_tot = size(pool_unmasked_corr, 2);
nCalcio_tot    = size(pool_calcio_raw, 2);

fprintf('  ROI unmasked pooled  : %i\n', nUnmasked_tot);
fprintf('  ROI solo-cloro pooled: %i\n\n', nCalcio_tot);

% ---- FIGURA 1: Media ROI Unmasked (corretta vs raw) ---------------------

mean_corr = median(pool_unmasked_corr, 2, 'omitnan');
std_corr = arrayfun(@(t) mad(pool_unmasked_corr(t, ~isnan(pool_unmasked_corr(t,:))), 1), ...
           (1:size(pool_unmasked_corr,1))');
mean_raw  = median(pool_unmasked_raw,  2, 'omitnan');
std_raw = arrayfun(@(t) mad(pool_unmasked_raw(t, ~isnan(pool_unmasked_raw(t,:))), 1), ...
          (1:size(pool_unmasked_raw,1))');

fig1_title = sprintf('ROI Unmasked (n=%i) | Stim: %s', ...
    nUnmasked_tot, strjoin(stimToAnalyze,'+'));

figure('Name', fig1_title, 'NumberTitle','off','Color','w');
tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');

% Subplot 1: DeltaGG_920 corretto (post-avatar)
nexttile; hold on
fill([Time_ref; flipud(Time_ref)], ...
     [mean_corr+std_corr; flipud(mean_corr-std_corr)], ...
     col_unm, 'FaceAlpha',0.25, 'EdgeColor','none');
plot(Time_ref, mean_corr, 'Color', col_unm, 'LineWidth', 2.0);
xline(0,'k--'); yline(0,'k:');
xlabel('Time (s)'); ylabel('\DeltaG/G (%)');
title(sprintf('\\DeltaG/G 920nm CORRETTO\nn=%i ROI unmasked', nUnmasked_tot));
set(gca,'Box','off');
xlim([Time_ref(1), Time_ref(end)]);

% Subplot 2: DeltaGG_920 raw (G_smooth non corretto)
nexttile; hold on
fill([Time_ref; flipud(Time_ref)], ...
     [mean_raw+std_raw; flipud(mean_raw-std_raw)], ...
     [0.4 0.4 0.8], 'FaceAlpha',0.25, 'EdgeColor','none');
plot(Time_ref, mean_raw, 'Color',[0.2 0.2 0.8], 'LineWidth', 2.0);
xline(0,'k--'); yline(0,'k:');
xlabel('Time (s)'); ylabel('\DeltaG/G (%)');
title(sprintf('\\DeltaG/G 920nm RAW (no correzione avatar)\nn=%i ROI unmasked', nUnmasked_tot));
set(gca,'Box','off');
xlim([Time_ref(1), Time_ref(end)]);

sgtitle(fig1_title, 'Interpreter','none');

% ---- FIGURA 2: Media ROI solo-cloro (raw) --------------------------------

mean_calcio = median(pool_calcio_raw, 2, 'omitnan');
std_calcio = arrayfun(@(t) mad(pool_calcio_raw(t, ~isnan(pool_calcio_raw(t,:))), 1), ...
             (1:size(pool_calcio_raw,1))');

fig2_title = sprintf('ROI Solo-Calcio (n=%i) | Stim: %s', ...
    nCalcio_tot, strjoin(stimToAnalyze,'+'));

figure('Name', fig2_title, 'NumberTitle','off','Color','w');
hold on
fill([Time_ref; flipud(Time_ref)], ...
     [mean_calcio+std_calcio; flipud(mean_calcio-std_calcio)], ...
     col_norm, 'FaceAlpha',0.20, 'EdgeColor','none');
plot(Time_ref, mean_calcio, 'Color', col_norm, 'LineWidth', 2.0);
xline(0,'k--'); yline(0,'k:');
xlabel('Time (s)'); ylabel('\DeltaG/G (%)');
title(sprintf('\\DeltaG/G 920nm RAW | Solo-Calcio (n=%i)', nCalcio_tot));
set(gca,'Box','off');
xlim([Time_ref(1), Time_ref(end)]);
sgtitle(fig2_title, 'Interpreter','none');

fprintf('Figure riassuntive generate.\n\n');
% ---- FIGURA 1b: Differenza (corretto - non corretto) --------------------
% Rappresenta il contributo dell'avatar sottratto dalla correzione.
% Se la correzione funziona bene, questa traccia dovrebbe avere la forma
% dell'avatar (valle negativa seguita da recovery) e ampiezza coerente
% tra le ROI unmasked.

% Le due pool hanno lo stesso numero di colonne (stesse ROI, stesso ordine)
% quindi la differenza è elemento per elemento.
diff_pool = pool_unmasked_corr - pool_unmasked_raw;  % nTime x nROIs

mean_diff = median(diff_pool, 2, 'omitnan');
std_diff  = arrayfun(@(t) mad(diff_pool(t, ~isnan(diff_pool(t,:))), 1), ...
            (1:size(diff_pool,1))');

fig1b_title = sprintf('Contributo Avatar (corretto - raw) | ROI Unmasked (n=%i) | Stim: %s', ...
    nUnmasked_tot, strjoin(stimToAnalyze,'+'));

figure('Name', fig1b_title, 'NumberTitle','off','Color','w');
hold on

fill([Time_ref; flipud(Time_ref)], ...
     [mean_diff+std_diff; flipud(mean_diff-std_diff)], ...
     [0.8 0.4 0.0], 'FaceAlpha',0.25, 'EdgeColor','none');
plot(Time_ref, mean_diff, 'Color',[0.7 0.3 0.0], 'LineWidth', 2.0);
yline(0,'k:');
xline(0,'k--');

xlabel('Time (s)');
ylabel('\DeltaG/G (%)');
title(sprintf('Contributo avatar (corretto - raw)\nn=%i ROI unmasked', nUnmasked_tot));
set(gca,'Box','off');
xlim([Time_ref(1), Time_ref(end)]);

sgtitle(fig1b_title, 'Interpreter','none');

fprintf('  Figura contributo avatar generata.\n\n');

% ---- FIGURA 1c: Media DeltaG/R (830 vs 920 nm) --------------------------
% Calcolo statistiche (Mediana e MAD) per G/R
mean_GR830 = median(pool_unmasked_GR830, 2, 'omitnan');
std_GR830  = arrayfun(@(t) mad(pool_unmasked_GR830(t, ~isnan(pool_unmasked_GR830(t,:))), 1), ...
             (1:size(pool_unmasked_GR830,1))');

mean_GR920 = median(pool_unmasked_GR920, 2, 'omitnan');
std_GR920  = arrayfun(@(t) mad(pool_unmasked_GR920(t, ~isnan(pool_unmasked_GR920(t,:))), 1), ...
             (1:size(pool_unmasked_GR920,1))');

fig1c_title = sprintf('Rapporto G/R Unmasked (n=%i) | Stim: %s', ...
    nUnmasked_tot, strjoin(stimToAnalyze,'+'));
figure('Name', fig1c_title, 'NumberTitle','off','Color','w');
tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');

% Subplot 1: DeltaG/R a 830nm
nexttile; hold on
fill([Time_ref; flipud(Time_ref)], ...
     [mean_GR830+std_GR830; flipud(mean_GR830-std_GR830)], ...
     [0.8 0.4 0.8], 'FaceAlpha',0.25, 'EdgeColor','none'); % Viola/Magenta
plot(Time_ref, mean_GR830, 'Color', [0.6 0.2 0.6], 'LineWidth', 2.0);
xline(0,'k--'); yline(0,'k:');
xlabel('Time (s)'); ylabel('\DeltaG/R (%)');
title('\DeltaG/R 830nm (iClima) RAW');
set(gca,'Box','off'); xlim([Time_ref(1), Time_ref(end)]);

% Subplot 2: DeltaG/R a 920nm
nexttile; hold on
fill([Time_ref; flipud(Time_ref)], ...
     [mean_GR920+std_GR920; flipud(mean_GR920-std_GR920)], ...
     [0.4 0.7 1.0], 'FaceAlpha',0.25, 'EdgeColor','none'); % Azzurro
plot(Time_ref, mean_GR920, 'Color', [0.2 0.5 0.8], 'LineWidth', 2.0);
xline(0,'k--'); yline(0,'k:');
xlabel('Time (s)'); ylabel('\DeltaG/R (%)');
title('\DeltaG/R 920nm (GCaMP + iClima) RAW');
set(gca,'Box','off'); xlim([Time_ref(1), Time_ref(end)]);

sgtitle(fig1c_title, 'Interpreter','none');
fprintf('  Figura DeltaG/R generata.\n\n');

% =========================================================================
%% SEZIONE 5 — SALVATAGGIO
% =========================================================================

risposta_salva = questdlg('Vuoi salvare ClCa aggiornata (con campo IsUnmasked)?', ...
    'Salvataggio', 'Sì', 'No', 'Sì');
if strcmp(risposta_salva, 'Sì')
    defaultName = sprintf('ClCa_%s_session5_%s.mat', ...
        strjoin(stimToAnalyze,'+'), datestr(now,'yyyymmdd_HHMM')); %#ok<TNOW1,DATST>
    [saveName, saveDir] = uiputfile('*.mat','Salva ClCa', defaultName);
    if ~isequal(saveName, 0)
        save(fullfile(saveDir, saveName), 'ClCa', '-v7.3');
        fprintf('ClCa salvata in:\n  %s\n', fullfile(saveDir, saveName));
    else
        fprintf('Salvataggio annullato.\n');
    end
else
    fprintf('ClCa disponibile nel workspace.\n');
end

fprintf('\n=========================================================\n');
fprintf('               Sessione 5 completata.                    \n');
fprintf('=========================================================\n');





%% FUNZIONI LOCALI
% =========================================================================

% -------------------------------------------------------------------------
function [stimLabel, stimFound] = assignStimulus(fname)
% ASSIGNSTIMULUS  Identifica il tipo di stimolo dal nome del file.
%
% Gerarchia di controllo (importante: CycleGrid PRIMA di Grid):
%   1. CycleGrid / DoubleGrid  -> 'CycleGrid'
%   2. Grid con orientamento   -> 'Grid_XXX' (es. 'Grid_0', 'Grid_90')
%   3. Nessuna corrispondenza  -> stimFound = false
%
% Gli orientamenti riconosciuti per Grid sono:
%   0, 45, 90, 135, 180, 270, 315
%   Per Grid315 è accettata anche la forma 'GRID-315'.
%
% Input:
%   fname     - nome del file (stringa)
%
% Output:
%   stimLabel - stringa con il nome dello stimolo ('CycleGrid', 'Grid_0', ...)
%   stimFound - true se trovata corrispondenza, false altrimenti
 
    stimLabel = '';
    stimFound = false;
 
    % ------------------------------------------------------------------
    % 1. CYCLEGRID
    %    Pattern riconosciuti (case-insensitive):
    %      CycleGrid0_90, CycleGrid, DoubleGrid0_90, DoubleGrid
    %    Nota: va controllato PRIMA di Grid perché il nome contiene
    %    anche "Grid", che verrebbe erroneamente matchato dalla regola 2.
    % ------------------------------------------------------------------
    if ~isempty(regexpi(fname, '(CycleGrid|DoubleGrid)', 'once'))
        stimLabel = 'CycleGrid';
        stimFound = true;
        return
    end
 
    % ------------------------------------------------------------------
    % 2. GRID con orientamento specifico
    %    Orientamenti riconosciuti: 0, 45, 90, 135, 180, 270, 315
    %
    %    Pattern riconosciuti per ciascun orientamento (es. per 90):
    %      Grid90, GRID90, G90, G90IsoOn, grid90, Grid-90
    %      G0_830nm, G90_130um, GRID-315, ecc.
    %
    %    NOTA sul perché NON usiamo \b (word boundary):
    %    In MATLAB, \b considera '_' come carattere word (fa parte di \w),
    %    quindi tra '0' e '_' NON c'è word boundary. Per esempio
    %    'G0_830nm' non viene matchato da 'G0\b'.
    %    Soluzione: lookahead (?=[_\-\.nm ]|$) che accetta esplicitamente
    %    i separatori reali nei nomi file: _ - . spazio nm oppure fine stringa.
    %    Questo cattura correttamente G0_830nm, G90_130um, G0IsoOn, ecc.
    % ------------------------------------------------------------------
    orientations = [0, 45, 90, 135, 180, 270, 315];
 
    for ang = orientations
        % Pattern: (Grid|G) + trattino opzionale + numero + lookahead separatore
        % Il lookahead (?=...) non consuma caratteri, serve solo a verificare
        % cosa segue il numero per evitare falsi match (es. G04 != G0).
        pattern = sprintf('(Grid|G)-?%d(?=[_\\-\\. ]|nm|IsoOn|$)', ang);
 
        if ~isempty(regexpi(fname, pattern, 'once'))
            stimLabel = sprintf('Grid_%d', ang);
            stimFound = true;
            return
        end
    end
 
    % ------------------------------------------------------------------
    % 3. Nessuna corrispondenza trovata
    % ------------------------------------------------------------------
    % stimLabel rimane '' e stimFound rimane false
 
end


% -------------------------------------------------------------------------
function smoothed = causalSmooth(data, win_frames)
% CAUSALSMOOTH  Applica una moving average causale ai dati.
%
% "Causale" significa che ogni punto viene mediato con i 'win_frames-1'
% punti precedenti (incluso il punto stesso), MAI con i punti futuri.
% Questo garantisce che lo smooth non introduca anticipazioni apparenti
% della risposta rispetto allo stimolo.
%
% Implementazione: usa movmean() di MATLAB con endpoint [win_frames-1, 0].
%   - win_frames-1 : quanti frame nel passato includere
%   - 0            : quanti frame nel futuro includere (zero = causale)
%
% I dati possono essere 1D, 2D o 3D (come G_sub_dark: frames x trials x rois).
% Lo smooth viene applicato lungo la prima dimensione (il tempo).
%
% Input:
%   data       - array numerico (prima dimensione = tempo)
%   win_frames - numero di frames della finestra (intero >= 1)
%
% Output:
%   smoothed   - array delle stesse dimensioni di data, smoothato

    if win_frames <= 1
        % Finestra di 1 frame = nessuno smooth (ogni punto è solo se stesso)
        smoothed = data;
        return
    end

    % movmean con endpoint [passato, futuro]:
    %   [win_frames-1, 0] -> media sui win_frames-1 frame precedenti
    %                        + il frame corrente, nessun frame futuro
    smoothed = movmean(data, [win_frames-1, 0], 1, 'omitnan');

end

%% FUNZIONE LOCALE
% =========================================================================

function idx_sorted = argsort_desc(values)
% Restituisce gli indici che ordinano 'values' in modo decrescente,
% mettendo i NaN alla fine.
    [~, idx_sorted] = sort(values, 'descend', 'MissingPlacement', 'last');
end