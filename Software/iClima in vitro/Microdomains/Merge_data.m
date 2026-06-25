%% analyze_microID.m
%Script for the analysis of the Cross-correlation among different
%stimulation trials, and for the analysis of the nearest neighbour for the
%chloride microdomain

%%
% Carica i risultati salvati da microID.mlx (file .mat) da cartelle Rec##,
% li organizza in una struct, e fa due plot di dataCD vs lastSh.
%
% Struttura cartelle attesa:
%   baseDir/
%     Rec01/   -> contiene uno o più file .mat con variabili dataCD, lastSh, nEvents, ...
%     Rec02/
%     ...
%
% Ogni .mat = una acquisizione (trial) di quel dendrite.

clear; clc; close all;

%% ── 1. PARAMETRI ─────────────────────────────────────────────────────────────
baseDir      = 'C:\Users\giak7\OneDrive\Documenti\LAB\New Chloride Sensor\Microdomini\Dataset_Thr3.5';      % cambia con il path della tua cartella base
edges        = 0:50;     % edges usati in microID per la CDF (0..20 pixel)
xAx          = edges(1:end);
nEventsRange = 61:65;    % indici di nEvents su cui sommare per la selezione
nEventsThr   = 5;        % soglia: acq scartata se sum(nEvents(61:65)) < nEventsThr

%% ── 2. LOOP SUI REC E CARICAMENTO ────────────────────────────────────────────
recDirs = dir(fullfile(baseDir, 'Rec*'));
recDirs = recDirs([recDirs.isdir]);

if isempty(recDirs)
    error('Nessuna cartella Rec* trovata in: %s', baseDir);
end

data = struct();

for r = 1:numel(recDirs)
    recName = recDirs(r).name;
    recPath = fullfile(baseDir, recName);

    matFiles = dir(fullfile(recPath, '*.mat'));
    if isempty(matFiles)
        warning('Nessun .mat in %s, salto.', recName);
        continue
    end

    data(r).recName = recName;
    data(r).acq     = struct();

    for t = 1:numel(matFiles)
        matPath = fullfile(recPath, matFiles(t).name);
        S = load(matPath);

        if ~isfield(S, 'dataCD') || ~isfield(S, 'lastSh')
            warning('File %s non contiene dataCD e/o lastSh, salto.', matFiles(t).name);
            continue
        end

        % Copia TUTTI i campi del .mat dentro data(r).acq(t)
        fields = fieldnames(S);
        data(r).acq(t).fileName = matFiles(t).name;
        for f = 1:numel(fields)
            data(r).acq(t).(fields{f}) = S.(fields{f});
        end
        % Assicura che dataCD e lastSh siano righe 1xN
        data(r).acq(t).dataCD  = S.dataCD(:)';
        data(r).acq(t).lastSh  = S.lastSh(:)';
        % Campo di selezione: inizialmente tutte incluse
        data(r).acq(t).include = true;
    end

    fprintf('Caricato %s: %d acquisizioni\n', recName, numel(data(r).acq));
end

fprintf('\nTotale dendrite caricati: %d\n', numel(data));

%% ── 3. DIALOGO thrZ + SELEZIONE ACQUISIZIONI ─────────────────────────────────

% Leggi thrZ da una qualsiasi acquisizione valida
thrZ = NaN;
for r = 1:numel(data)
    for t = 1:numel(data(r).acq)
        if isfield(data(r).acq(t), 'thrZ') && ~isempty(data(r).acq(t).thrZ)
            thrZ = data(r).acq(t).thrZ;
            break
        end
    end
    if ~isnan(thrZ), break, end
end

% Finestra di dialogo
if isnan(thrZ)
    msg = sprintf(['thrZ non trovato nei dati.\n\n' ...
                   'Vuoi applicare il criterio di selezione basato su nEvents?\n' ...
                   '(somma nEvents(61:65) >= %d)'], nEventsThr);
else
    msg = sprintf(['thrZ letto dai dati: %.4f\n\n' ...
                   'Vuoi applicare il criterio di selezione basato su nEvents?\n' ...
                   '(somma nEvents(61:65) >= %d)'], thrZ, nEventsThr);
end

answer = questdlg(msg, 'Selezione acquisizioni', 'Sì', 'No', 'Sì');
useSelection = strcmp(answer, 'Sì');

% Applica selezione (o non)
nScartate = 0;
if useSelection
    fprintf('\n── Criterio di selezione ATTIVO (soglia nEvents = %d) ──\n', nEventsThr);
    for r = 1:numel(data)
        for t = 1:numel(data(r).acq)
            acq = data(r).acq(t);
            if ~isfield(acq, 'nEvents') || isempty(acq.nEvents)
                warning('%s / %s: nEvents mancante, acq esclusa.', ...
                        data(r).recName, acq.fileName);
                data(r).acq(t).include = false;
                nScartate = nScartate + 1;
                continue
            end
            nEv = acq.nEvents(:);
            % Clamp degli indici al range disponibile
            idx = nEventsRange(nEventsRange <= numel(nEv));
            if isempty(idx) || sum(nEv(idx)) < nEventsThr
                data(r).acq(t).include = false;
                nScartate = nScartate + 1;
                fprintf('  SCARTATA: %s / %s  (sum nEvents(61:65) = %d)\n', ...
                        data(r).recName, acq.fileName, sum(nEv(idx)));
            end
        end
    end
    fprintf('Acquisizioni scartate: %d\n', nScartate);
else
    fprintf('\n── Criterio di selezione NON applicato: tutte le acq incluse ──\n');
end

%% ── 4. RACCOLTA GLOBALE PER I PLOT ───────────────────────────────────────────
% Solo le acquisizioni con include=true entrano nelle matrici e nelle medie.

allDataCD  = [];
allLastSh  = [];
trialRec   = [];

meanDataCD_perRec = zeros(numel(data), numel(xAx));
meanLastSh_perRec = zeros(numel(data), numel(xAx));

for r = 1:numel(data)
    nTrials = numel(data(r).acq);
    recDataCD = [];
    recLastSh = [];

    for t = 1:nTrials
        if ~data(r).acq(t).include, continue, end

        cd_t = data(r).acq(t).dataCD;
        sh_t = data(r).acq(t).lastSh;
        n    = min(numel(cd_t), numel(xAx));

        row_cd = NaN(1, numel(xAx)); row_cd(1:n) = cd_t(1:n);
        row_sh = NaN(1, numel(xAx)); row_sh(1:n) = sh_t(1:n);

        recDataCD = [recDataCD; row_cd]; %#ok<AGROW>
        recLastSh = [recLastSh; row_sh]; %#ok<AGROW>
    end

    if isempty(recDataCD)
        meanDataCD_perRec(r,:) = NaN;
        meanLastSh_perRec(r,:) = NaN;
        fprintf('ATTENZIONE: %s non ha acquisizioni valide dopo la selezione.\n', data(r).recName);
        continue
    end

    allDataCD = [allDataCD; recDataCD]; %#ok<AGROW>
    allLastSh = [allLastSh; recLastSh]; %#ok<AGROW>
    trialRec  = [trialRec; r * ones(size(recDataCD,1),1)]; %#ok<AGROW>

    meanDataCD_perRec(r,:) = mean(recDataCD, 1, 'omitnan');
    meanLastSh_perRec(r,:) = mean(recLastSh, 1, 'omitnan');
end

nIncluded = size(allDataCD, 1);
fprintf('\nAcquisizioni incluse nell''analisi: %d\n', nIncluded);

%% ── 5. COLORI ────────────────────────────────────────────────────────────────
nRec   = numel(data);
cmap   = lines(nRec);

boldColorCD = [0.12 0.47 0.71];   % blu
boldColorSh = [0.84 0.15 0.16];   % rosso
alphaLine   = 0.25;
alphaRec    = 0.35;

%% ── 6. PLOT 1 – TUTTE LE ACQUISIZIONI DI TUTTI I REC ────────────────────────
fig1 = figure('Name','Plot 1 – Tutti i trial, tutti i Rec', ...
              'Position',[100 100 900 500]);
hold on; box on; grid on;

for i = 1:size(allDataCD, 1)
    r = trialRec(i);
    c = cmap(r,:);
    plot(xAx, allDataCD(i,:), '-',  'Color', [c alphaLine], 'LineWidth',0.8);
    plot(xAx, allLastSh(i,:), '--', 'Color', [c alphaLine], 'LineWidth',0.8);
end

globalMeanCD = mean(allDataCD, 1, 'omitnan');
globalMeanSh = mean(allLastSh, 1, 'omitnan');
plot(xAx, globalMeanCD, '-',  'Color', boldColorCD, 'LineWidth',3);
plot(xAx, globalMeanSh, '--', 'Color', boldColorSh, 'LineWidth',3);

xlabel('Nearest-neighbour distance (pixel)');
ylabel('Cumulative probability');
selStr = 'selezione ON'; if ~useSelection, selStr = 'selezione OFF'; end
title(sprintf('Plot 1 – Tutti i trial, tutti i dendrite  [%s, n=%d]', selStr, nIncluded));

h_dummy_cd = plot(NaN, NaN, '-',  'Color', boldColorCD, 'LineWidth',2);
h_dummy_sh = plot(NaN, NaN, '--', 'Color', boldColorSh, 'LineWidth',2);
legend([h_dummy_cd, h_dummy_sh], {'dataCD (media globale)', 'lastSh (media globale)'}, ...
       'Location','southeast');
ylim([0 1]); xlim([xAx(1) xAx(end)]);

%% ── 7. PLOT 2 – UNA LINEA PER REC (media delle acq) ─────────────────────────
fig2 = figure('Name','Plot 2 – Media per Rec', ...
              'Position',[150 150 900 500]);
hold on; box on; grid on;

for r = 1:nRec
    plot(xAx, meanDataCD_perRec(r,:), '-',  'Color', [boldColorCD alphaRec], 'LineWidth',1.5);
    plot(xAx, meanLastSh_perRec(r,:), '--', 'Color', [boldColorSh alphaRec], 'LineWidth',1.5);
end

grandMeanCD = mean(meanDataCD_perRec, 1, 'omitnan');
grandMeanSh = mean(meanLastSh_perRec, 1, 'omitnan');
plot(xAx, grandMeanCD, '-',  'Color', boldColorCD, 'LineWidth',3);
plot(xAx, grandMeanSh, '--', 'Color', boldColorSh, 'LineWidth',3);

xlabel('Nearest-neighbour distance (pixel)');
ylabel('Cumulative probability');
title(sprintf('Plot 2 – Media per dendrite  [%s]', selStr));

h_dummy_cd = plot(NaN, NaN, '-',  'Color', boldColorCD, 'LineWidth',2.5);
h_dummy_sh = plot(NaN, NaN, '--', 'Color', boldColorSh, 'LineWidth',2.5);
legend([h_dummy_cd, h_dummy_sh], {'dataCD (grand mean)', 'lastSh (grand mean)'}, ...
       'Location','southeast');
ylim([0 1]); xlim([xAx(1) xAx(end)]);

%% ── 8. STATISTICA ────────────────────────────────────────────────────────────
distPx  = 5;
distIdx = distPx + 1;

fprintf('\n══════════════════════════════════════════════════\n');
fprintf('  STATISTICA: dataCD vs lastSh  [%s]\n', selStr);
fprintf('══════════════════════════════════════════════════\n');

scalarCD_A = allDataCD(:, distIdx);
scalarSh_A = allLastSh(:, distIdx);
scalarCD_B = trapz(xAx, allDataCD, 2);
scalarSh_B = trapz(xAx, allLastSh, 2);

% Livello TRIAL
fprintf('\n── Livello TRIAL (n = %d) ──\n', nIncluded);
[p_A, h_A, stats_A] = signrank(scalarCD_A, scalarSh_A);
fprintf('  Opzione A (CDF a dist=%d px):  mediana diff = %+.4f,  p = %.4f,  z = %.3f,  H = %d\n', ...
        distPx, median(scalarCD_A - scalarSh_A,'omitnan'), p_A, stats_A.zval, h_A);

[p_B, h_B, stats_B] = signrank(scalarCD_B, scalarSh_B);
fprintf('  Opzione B (AUC):               mediana diff = %+.4f,  p = %.4f,  z = %.3f,  H = %d\n', ...
        median(scalarCD_B - scalarSh_B,'omitnan'), p_B, stats_B.zval, h_B);

% Livello REC
recMeanCD_A = NaN(nRec,1); recMeanSh_A = NaN(nRec,1);
recMeanCD_B = NaN(nRec,1); recMeanSh_B = NaN(nRec,1);
for r = 1:nRec
    idx = trialRec == r;
    if any(idx)
        recMeanCD_A(r) = mean(scalarCD_A(idx),'omitnan');
        recMeanSh_A(r) = mean(scalarSh_A(idx),'omitnan');
        recMeanCD_B(r) = mean(scalarCD_B(idx),'omitnan');
        recMeanSh_B(r) = mean(scalarSh_B(idx),'omitnan');
    end
end

fprintf('\n── Livello REC (n = %d dendrite) ──\n', nRec);
if nRec >= 6
    [p_rA, h_rA, s_rA] = signrank(recMeanCD_A, recMeanSh_A);
    fprintf('  Opzione A (CDF a dist=%d px):  mediana diff = %+.4f,  p = %.4f,  z = %.3f,  H = %d\n', ...
            distPx, median(recMeanCD_A - recMeanSh_A,'omitnan'), p_rA, s_rA.zval, h_rA);
    [p_rB, h_rB, s_rB] = signrank(recMeanCD_B, recMeanSh_B);
    fprintf('  Opzione B (AUC):               mediana diff = %+.4f,  p = %.4f,  z = %.3f,  H = %d\n', ...
            median(recMeanCD_B - recMeanSh_B,'omitnan'), p_rB, s_rB.zval, h_rB);
else
    fprintf('  n = %d Rec (n basso, test indicativo).\n', nRec);
    fprintf('  Opzione A:  mediana diff = %+.4f\n', median(recMeanCD_A - recMeanSh_A,'omitnan'));
    fprintf('  Opzione B:  mediana diff = %+.4f\n', median(recMeanCD_B - recMeanSh_B,'omitnan'));
end

fprintf('\n══════════════════════════════════════════════════\n');
fprintf('Fatto.\n');


%% ── 9. PLOT DISTRIBUZIONE FWHM (Width) ────────────────────────────────────────
allWidths = []; % Inizializza il contenitore per il pooling

for r = 1:numel(data)
    for t = 1:numel(data(r).acq)
        % Consideriamo solo le acquisizioni incluse nella selezione
        if data(r).acq(t).include
            % Controlliamo se la struttura 'domain' e il campo 'width' esistono
            if isfield(data(r).acq(t), 'domain') && isfield(data(r).acq(t).domain, 'width')
                w = data(r).acq(t).domain.width;
                if ~isempty(w)
                    % Aggiungiamo i valori (tutti, se ce n'è più di uno)
                    allWidths = [allWidths; w(:)]; %#ok<AGROW>
                end
            end
        end
    end
end

if isempty(allWidths)
    warning('Nessun valore di FWHM (domain.width) trovato nelle acquisizioni incluse.');
else
    % Calcolo statistiche descrittive
    medW = median(allWidths);
    avgW = mean(allWidths);
    stdW = std(allWidths);

    fig3 = figure('Name','Distribuzione FWHM Microdomini', 'Position',[200 200 600 500]);
    hold on; box on; grid on;
    
    % Istogramma con densità (Normalization 'pdf' o 'count' a tua scelta)
    h = histogram(allWidths,150, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'w');
   
    
    % Aggiunta linee per media e mediana
    yLimits = ylim;
    line([medW medW], yLimits, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    line([avgW avgW], yLimits, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '-.');
    
    xlabel('FWHM - Full Width at Half Maximum');
    ylabel('Conteggio (n)');
    title(sprintf('Distribuzione Spaziale Microdomini (n = %d)', numel(allWidths)));
    
    legend('FWHM osservati', ...
           ['Mediana: ' num2str(medW, '%.2f')], ...
           ['Media: ' num2str(avgW, '%.2f')], ...
           'Location', 'northeast');
    
    fprintf('\n── ANALISI LARGHEZZA MICRODOMINI ──\n');
    fprintf('  Totale domini analizzati: %d\n', numel(allWidths));
    fprintf('  Media FWHM:   %.3f\n', avgW);
    fprintf('  Mediana FWHM: %.3f\n', medW);
    fprintf('  Dev. Std:     %.3f\n', stdW);
end

%% ── p_val per trial (da microID) ────────────────────────────────────────────
% Raccoglie tutti i p_val salvati in data.acq.p_val e conta quanti
% sono significativi (< 0.05) sul totale dei trial inclusi.

pValThr = 0.05;

allPval   = [];
recLabels_pval = {};

for r = 1:numel(data)
    for t = 1:numel(data(r).acq)
        if ~data(r).acq(t).include, continue, end
        if ~isfield(data(r).acq(t), 'p_val') || isempty(data(r).acq(t).p_val)
            warning('%s / %s: p_val mancante, salto.', data(r).recName, data(r).acq(t).fileName);
            continue
        end
        allPval(end+1) = data(r).acq(t).p_val; %#ok<AGROW>
        recLabels_pval{end+1} = data(r).recName; %#ok<AGROW>
    end
end

nTotal = numel(allPval);
nSig   = sum(allPval < pValThr);

fprintf('\n── p_val da microID (threshold = %.2f) ──\n', pValThr);
fprintf('  Trial inclusi con p_val:  %d\n',   nTotal);
fprintf('  Significativi:            %d / %d  (%.1f%%)\n', nSig, nTotal, 100*nSig/nTotal);

% Breakdown per Rec
fprintf('\n  Dettaglio per dendrite:\n');
for r = 1:numel(data)
    idx = strcmp(recLabels_pval, data(r).recName);
    if ~any(idx), continue, end
    pv_r   = allPval(idx);
    nS_r   = sum(pv_r < pValThr);
    fprintf('    %s:  %d / %d  significativi  (%.1f%%)\n', ...
            data(r).recName, nS_r, numel(pv_r), 100*nS_r/numel(pv_r));
end


%% ── FOREST PLOT p_val per acquisizione ──────────────────────────────────────
sigThr   = 0.05;
log10Thr = -log10(sigThr);   % 1.301

% Numero massimo di acquisizioni tra tutti i Rec
nAcqMax = 0;
for r = 1:numel(data)
    nAcqMax = max(nAcqMax, numel(data(r).acq));
end

% Matrice pVal: righe = Rec (1=Rec01, fine=RecN), colonne = Acq
pValMat = NaN(numel(data), nAcqMax);
inclMat = false(numel(data), nAcqMax);

for r = 1:numel(data)
    for t = 1:numel(data(r).acq)
        % Ricava l'indice della colonna dal nome file (Acq1→1, Acq2→2, Acq3→3)
        fname = data(r).acq(t).fileName;
        acqNum = str2double(regexp(fname, '(?<=Acq)\d+', 'match', 'once'));
        if isnan(acqNum), continue, end

        pValMat(r, acqNum) = data(r).acq(t).p_val;
        inclMat(r, acqNum) = data(r).acq(t).include;
    end
end
logPmat = -log10(pValMat);

acqColors = [
    0.18 0.45 0.69;
    0.90 0.55 0.00;
    0.76 0.15 0.15;
];
if nAcqMax > 3
    acqColors = [acqColors; lines(nAcqMax - 3)];
end

nRec_    = numel(data);
recNames = {data.recName};

% yPos crescente (1=Rec01 in basso, nRec=RecN in cima) — MATLAB vuole valori crescenti
yPos     = 1:nRec_;
% Etichette: Rec01 in basso → l'ordine è già 1..N
yLabels  = recNames;   % recNames{1}='Rec01' va a y=1 (basso), recNames{end} a y=nRec_ (cima)

xMax = max(logPmat(:), [], 'omitnan');
xMax = max(xMax, log10Thr + 0.5);

figP = figure('Name','Forest plot – p_val per acquisizione', ...
              'Position',[200 50 300*nAcqMax 80 + 20*nRec_]);

for acq = 1:nAcqMax
    ax = subplot(1, nAcqMax, acq);
    hold on; box on;

    col = acqColors(min(acq, size(acqColors,1)), :);

    for r = 1:nRec_
        y = yPos(r);

        if isnan(pValMat(r, acq))
            % Acq mancante: piccolo segno grigio
            plot(0, y, '-', 'Color', [0.88 0.88 0.88], 'LineWidth', 0.5);
            continue
        end

        lp    = logPmat(r, acq);
        isSig = pValMat(r, acq) < sigThr;

        if isSig
            dotCol = col;  mkSize = 7;
        else
            dotCol = [0.65 0.65 0.65];  mkSize = 5;
        end

        % Punto
        plot(lp, y, 'd', 'MarkerFaceColor', dotCol, ...
             'MarkerEdgeColor', dotCol, 'MarkerSize', mkSize);

        % Numero del p_val accanto al punto
        text(lp + 0.05, y, sprintf('%.3f', pValMat(r,acq)), ...
             'FontSize', 6, 'VerticalAlignment', 'middle', 'Color', dotCol);
    end

    % Linea soglia
    xline(log10Thr, '--k', 'p=0.05', ...
          'LabelVerticalAlignment','bottom', ...
          'LabelHorizontalAlignment','left', ...
          'FontSize', 7);

    % Asse Y: etichette Rec, ordine crescente (Rec01 in basso)
    set(ax, 'YTick', yPos, 'YTickLabel', yLabels, ...
            'YLim', [0.5, nRec_+0.5], ...
            'FontSize', 8);

    if acq > 1
        set(ax, 'YTickLabel', []);
    end

    xlabel('-log_{10}(p)', 'FontSize', 9);
    title(sprintf('Acq%d', acq), 'FontSize', 10);
    xlim([0, xMax]);
    grid on;
end

sgtitle('p\_val per dendrite e acquisizione  (-log_{10} scale)', 'FontSize', 11);