%Generation of a structure with all the iClima/GCaMP6f data
%  =========================================================================
% build_ClCa.m
% =========================================================================
%
% DESCRIZIONE:
%   Questo script costruisce (o aggiorna) la struttura dati "ClCa", che
%   raccoglie e organizza in coppie i dati di imaging a due lunghezze
%   d'onda (830 nm e 920 nm) acquisiti dallo stesso replicato biologico.
%
% STRUTTURA DI ClCa:
%   Ogni riga di ClCa rappresenta una coppia di esperimenti (830nm + 920nm)
%   dello stesso campo e dello stesso animale. I campi principali sono:
%       - Date, Mouse, Field, Stimulus (vuoto, da riempire in seguito)
%       - FoldName        : cartella di provenienza
%       - FileName_830    : nome del file a 830 nm
%       - FileName_920    : nome del file a 920 nm
%       - G_sub_dark_830/920   : canale verde con dark sottratto
%       - R_sub_dark_830/920   : canale rosso con dark sottratto
%       - Time_830/920         : vettore dei tempi
%       - Stimulus_time_830/920: tempi degli stimoli
%
% UTILIZZO:
%   1. Avvia lo script da MATLAB.
%   2. Ti verrà chiesto se vuoi partire da una struttura ClCa esistente
%      (per aggiungere nuovi dati) o crearne una da zero.
%   3. A ogni ciclo, si apre una finestra di dialogo: seleziona tutti i
%      file *_analyzed.mat di UNA cartella (quelli della stessa sessione).
%      Puoi selezionare più file contemporaneamente con Ctrl+click.
%   4. Premi "Annulla" (Cancel) nella finestra di dialogo per terminare
%      il ciclo di caricamento e salvare la struttura finale.
%   5. Lo script accoppia automaticamente i file 830 nm e 920 nm,
%      estrae i metadati dal nome della cartella e salva ClCa.
%
% FORMATO ATTESO DEL NOME CARTELLA:
%   YYYYMMDD_MouseXX_HHMM
%   L'identificativo del mouse può essere numerico, alfabetico o
%   alfanumerico (es. Mouse17, MouseD, MouseA3).
%
% FORMATO ATTESO DEI NOMI FILE:
%   I file devono contenere "830" o "920" nel nome.
%   Differenze tollerate tra i due file di una coppia:
%     - lunghezza d'onda (830nm vs 920nm)
%     - potenza (es. 30mW vs 20mW)
%     - numero progressivo di acquisizione (es. -004- vs -001-)
%
% NOTE:
%   - I file devono contenere una variabile MATLAB chiamata "out_struct".
%   - Se un file non può essere accoppiato, viene segnalato in console.
%   - Se il nome del file non contiene "C" maiuscola seguita da numero,
%     la coppia viene saltata e segnalata in console.
%
% =========================================================================

clearvars -except ClCa   % pulisce tutto tranne un eventuale ClCa pre-esistente
clc

fprintf('=========================================================\n');
fprintf('            BUILD_CLCA - Avvio dello script              \n');
fprintf('=========================================================\n\n');

% -------------------------------------------------------------------------
%% SEZIONE 1: Chiede se esiste già una struttura ClCa di partenza
% -------------------------------------------------------------------------

risposta = questdlg( ...
    'Esiste già una struttura ClCa da cui partire?', ...
    'Struttura esistente', ...
    'Sì, caricala', 'No, parti da zero', 'Sì, caricala');

switch risposta

    case 'Sì, caricala'
        % L'utente vuole partire da una ClCa già esistente
        [fname_clca, fdir_clca] = uigetfile('*.mat', ...
            'Seleziona il file .mat contenente ClCa');
        if isequal(fname_clca, 0)
            error('Nessun file selezionato. Script terminato.');
        end
        tmp = load(fullfile(fdir_clca, fname_clca), 'ClCa');
        if ~isfield(tmp, 'ClCa')
            error('Il file selezionato non contiene una variabile "ClCa".');
        end
        ClCa = tmp.ClCa;
        fprintf('Struttura ClCa caricata: %i righe esistenti.\n\n', numel(ClCa));

    case 'No, parti da zero'
        % Parte da una struttura vuota
        ClCa = struct();
        ClCa(:) = [];   % array vuoto con la struttura corretta
        fprintf('Partenza da struttura vuota.\n\n');

    otherwise
        % L'utente ha chiuso la finestra senza scegliere
        error('Nessuna scelta effettuata. Script terminato.');
end


% -------------------------------------------------------------------------
%% SEZIONE 2: Loop di caricamento file
% -------------------------------------------------------------------------
% Ad ogni iterazione l'utente seleziona UN BATCH di file dalla STESSA
% cartella. Premendo "Annulla" il loop si interrompe.

lastFolder = pwd;   % ricorda l'ultima cartella visitata
nRoundTot  = 0;     % contatore dei round completati

fprintf('--- Inizio ciclo di selezione file ---\n');
fprintf('Seleziona i file *_analyzed.mat dalla finestra di dialogo.\n');
fprintf('Premi ANNULLA per terminare il ciclo.\n\n');

while true

    % ---- Apre la finestra di dialogo per la selezione dei file ----------
    % Mostra solo i file che terminano con "_analyzed.mat"
    [fnames, fdir] = uigetfile( ...
        fullfile(lastFolder, '*_analyzed.mat'), ...
        'Seleziona file _analyzed.mat (Annulla = fine ciclo)', ...
        'MultiSelect', 'on');

    % Se l'utente preme Annulla, uigetfile restituisce 0
    if isequal(fnames, 0)
        fprintf('\nAnnulla premuto: fine del ciclo di caricamento.\n\n');
        break
    end

    % uigetfile restituisce una cella se si selezionano più file,
    % una stringa se ne viene selezionato uno solo: normalizziamo a cella
    if ischar(fnames)
        fnames = {fnames};
    end

    nFiles     = numel(fnames);
    nRoundTot  = nRoundTot + 1;
    lastFolder = fdir;   % aggiorna la cartella di riferimento

    fprintf('Round %i: selezionati %i file da:\n  %s\n', ...
        nRoundTot, nFiles, fdir);

    % ---- Estrae metadati dal nome della cartella ------------------------
    % Formato atteso: YYYYMMDD_MouseXX_HHMM
    % Prende la penultima cartella nel path (la cartella dell'esperimento)
    folderParts = strsplit(strtrim(fdir), filesep);
    folderParts(cellfun(@isempty, folderParts)) = []; % rimuove parti vuote
    sessionFolder = folderParts{end};  % nome della cartella sessione

    % Estrae data e mouse dal nome della cartella con espressione regolare
    % Pattern: 8 cifre per la data, underscore, "Mouse" + identificativo
    % L'identificativo del mouse può essere numerico (es. Mouse17) oppure
    % alfabetico (es. MouseD) o alfanumerico (es. MouseA3).
    tok_date  = regexp(sessionFolder, '^(\d{8})', 'tokens', 'once');
    tok_mouse = regexp(sessionFolder, '[Mm]ouse([A-Za-z0-9]+)', 'tokens', 'once');

    if ~isempty(tok_date)
        sessionDate = tok_date{1};   % es. '20260420'
    else
        sessionDate = 'unknown';
        warning('Impossibile estrarre la data dalla cartella: %s', sessionFolder);
    end

    if ~isempty(tok_mouse)
        sessionMouse = ['Mouse', tok_mouse{1}];   % es. 'Mouse17'
    else
        sessionMouse = 'unknown';
        warning('Impossibile estrarre il numero del mouse dalla cartella: %s', sessionFolder);
    end

    fprintf('  Data: %s | Mouse: %s\n', sessionDate, sessionMouse);

    % ---- Separa i file in 830 nm e 920 nm --------------------------------
    % Un file appartiene a 830 nm se contiene '830' nel nome, a 920 nm
    % se contiene '920'. (Con o senza 'nm'.)
    idx_830 = find(cellfun(@(n) ~isempty(regexp(n, '830', 'once')), fnames));
    idx_920 = find(cellfun(@(n) ~isempty(regexp(n, '920', 'once')), fnames));

    files_830 = fnames(idx_830);
    files_920 = fnames(idx_920);

    fprintf('  File a 830 nm: %i | File a 920 nm: %i\n', ...
        numel(files_830), numel(files_920));

    if isempty(files_830) || isempty(files_920)
        warning('Round %i: trovati file solo per una lunghezza d''onda. Nessuna coppia creata.', nRoundTot);
        continue
    end

    % ---- Accoppia i file 830-920 ----------------------------------------
    % Strategia: per ogni file 830, cerca il corrispondente 920 il cui
    % nome "normalizzato" (senza la potenza e senza la lunghezza d'onda)
    % è uguale.

    pairs = matchPairs(files_830, files_920);
    % 'pairs' è una matrice Nx2: colonna 1 = file 830, colonna 2 = file 920
    % Le righe con NaN indicano file non accoppiati.

    if isempty(pairs)
        warning('Round %i: nessuna coppia trovata.', nRoundTot);
        continue
    end

    fprintf('  Coppie trovate: %i\n', size(pairs, 1));

    % ---- Elabora ogni coppia e aggiunge a ClCa --------------------------
    for iPair = 1 : size(pairs, 1)

        fname830 = pairs{iPair, 1};
        fname920 = pairs{iPair, 2};

        fprintf('  Elaboro coppia %i/%i:\n    830: %s\n    920: %s\n', ...
            iPair, size(pairs,1), fname830, fname920);

        % -- Estrai il 'Field' dal nome del file 830 (vale anche per 920) --
        [fieldName, isValid] = extractField(fname830);

        if ~isValid
            % Il nome non rispetta il pattern C+numero: segnala e salta
            fprintf(['  [ATTENZIONE] Impossibile estrarre Field dal file:\n' ...
                '    830: %s\n    920: %s\n' ...
                '  La coppia viene SALTATA.\n'], fname830, fname920);
            continue
        end

        % -- Carica i dati dai file .mat -----------------------------------
        try
            data830 = load(fullfile(fdir, fname830), 'out_struct');
            s830    = data830.out_struct;
        catch ME
            warning('Errore nel caricare %s:\n%s', fname830, ME.message);
            continue
        end

        try
            data920 = load(fullfile(fdir, fname920), 'out_struct');
            s920    = data920.out_struct;
        catch ME
            warning('Errore nel caricare %s:\n%s', fname920, ME.message);
            continue
        end

        % -- Costruisce la nuova riga della struttura ClCa -----------------
        newRow = buildClCaRow( ...
            sessionDate, sessionMouse, fieldName, ...
            sessionFolder, fname830, fname920, ...
            s830, s920);

        % -- Aggiunge la riga a ClCa ---------------------------------------
        if isempty(ClCa) || (isstruct(ClCa) && numel(ClCa) == 0)
            ClCa = newRow;
        else
            ClCa(end+1) = newRow; %#ok<AGROW>
        end

    end % fine loop coppie

    fprintf('  ClCa ora contiene %i righe totali.\n\n', numel(ClCa));

end % fine while (loop round)


% -------------------------------------------------------------------------
%% SEZIONE 3: Salvataggio della struttura ClCa
% -------------------------------------------------------------------------

if numel(ClCa) == 0
    fprintf('Nessun dato caricato. ClCa è vuota, nessun file salvato.\n');
else
    fprintf('Totale righe in ClCa: %i\n', numel(ClCa));

    % Propone un nome file con data e ora correnti
    defaultName = ['ClCa_', datestr(now, 'yyyymmdd_HHMM'), '.mat']; %#ok<TNOW1,DATST>
    [saveName, saveDir] = uiputfile('*.mat', ...
        'Salva la struttura ClCa', defaultName);

    if isequal(saveName, 0)
        fprintf('Salvataggio annullato. ClCa è disponibile nel workspace.\n');
    else
        save(fullfile(saveDir, saveName), 'ClCa');
        fprintf('ClCa salvata in:\n  %s\n', fullfile(saveDir, saveName));
    end
end

fprintf('\n=========================================================\n');
fprintf('                   Script completato.                    \n');
fprintf('=========================================================\n');


% =========================================================================
%% FUNZIONI LOCALI
% =========================================================================

% -------------------------------------------------------------------------
function pairs = matchPairs(files_830, files_920)
% MATCHPAIRS  Accoppia file 830nm con file 920nm confrontando i nomi
%             normalizzati (senza lunghezza d'onda e senza potenza).
%
% Input:
%   files_830  - cell array di nomi file contenenti '830'
%   files_920  - cell array di nomi file contenenti '920'
%
% Output:
%   pairs      - cell array Nx2: {nome830, nome920} per ogni coppia trovata

    pairs = {};

    % Normalizza tutti i nomi: rimuove la lunghezza d'onda e la potenza
    norm_830 = cellfun(@normalizeName, files_830, 'UniformOutput', false);
    norm_920 = cellfun(@normalizeName, files_920, 'UniformOutput', false);

    usedIdx_920 = false(size(files_920)); % tiene traccia degli 920 già usati

    for i = 1 : numel(files_830)
        matched = false;
        for j = 1 : numel(files_920)
            if usedIdx_920(j)
                continue   % questo file 920 è già stato usato
            end
            if strcmp(norm_830{i}, norm_920{j})
                pairs(end+1, :) = {files_830{i}, files_920{j}}; %#ok<AGROW>
                usedIdx_920(j)  = true;
                matched = true;
                break
            end
        end
        if ~matched
            fprintf('  [ATTENZIONE] Nessuna coppia trovata per:\n    %s\n', files_830{i});
        end
    end

    % Segnala eventuali file 920 rimasti senza coppia
    for j = 1 : numel(files_920)
        if ~usedIdx_920(j)
            fprintf('  [ATTENZIONE] Nessuna coppia trovata per:\n    %s\n', files_920{j});
        end
    end
end


% -------------------------------------------------------------------------
function normName = normalizeName(fname)
% NORMALIZENAME  Rimuove dal nome file le parti variabili tra le due
%                lunghezze d'onda, in modo che i due file di una coppia
%                risultino identici dopo normalizzazione.
%
% Parti rimosse:
%   1. Estensione .mat
%   2. Lunghezza d'onda: '830nm', '920nm', '830', '920'  → 'WL'
%   3. Potenza: es. '30mW', '100mW'                      → 'PWR'
%   4. Numero progressivo di acquisizione: sequenza di 3 cifre preceduta
%      da '-' e seguita da '-' (es. '-004-', '-001-')    → '-NNN-'
%      Questo numero cambia tra i file 830 e 920 della stessa sessione
%      ma non è significativo per l'accoppiamento.
%
% Input:
%   fname    - stringa con il nome del file (con o senza estensione)
%
% Output:
%   normName - nome normalizzato, tutto in minuscolo

    % 1. Rimuove l'estensione .mat se presente
    normName = regexprep(fname, '\.mat$', '', 'ignorecase');

    % 2. Rimuove la lunghezza d'onda (prima con 'nm', poi senza,
    %    per evitare di lasciare residui come '-nm-')
    normName = regexprep(normName, '(830|920)nm', 'WL', 'ignorecase');
    normName = regexprep(normName, '(830|920)',    'WL', 'ignorecase');

    % 3. Rimuove la potenza: numero (1+ cifre) seguito da 'mW'
    normName = regexprep(normName, '\d+mW', 'PWR', 'ignorecase');

    % 4. Rimuove il numero progressivo di acquisizione.
    %    Questo numero è l'ultimo gruppo di cifre che precede il suffisso
    %    '_analyzed' (con eventuale trattino prima).
    %    Esempi: '-004-_analyzed' -> '-NNN-_analyzed'
    %            '-010-_analyzed' -> '-NNN-_analyzed'
    normName = regexprep(normName, '-\d+-_analyzed', '-NNN-_analyzed', 'ignorecase');

    % 5. Converte in minuscolo per confronto case-insensitive
    normName = lower(normName);
end


% -------------------------------------------------------------------------
function [fieldName, isValid] = extractField(fname)
% EXTRACTFIELD  Estrae il nome del campo (Field) dal nome del file.
%
% Regole:
%   - Cerca la prima 'C' maiuscola seguita da uno o più cifre.
%   - Dopo il numero deve esserci '-' o '_'. Se manca, viene aggiunto '_'.
%   - Se la 'C' non è seguita da cifre, isValid = false.
%
% Input:
%   fname     - stringa con il nome del file
%
% Output:
%   fieldName - stringa con il Field (es. 'C1', 'C3'), vuota se non valido
%   isValid   - true se il Field è stato estratto correttamente

    fieldName = '';
    isValid   = false;

    % Cerca "C" maiuscola seguita da cifre
    % Il token cattura la lettera C e le cifre che la seguono
    tok = regexp(fname, 'C(\d+)', 'tokens', 'once');

    if isempty(tok)
        % Nessuna "C" seguita da numero trovata
        return
    end

    fieldName = ['C', tok{1}];   % es. 'C1', 'C13'
    isValid   = true;

    % Controlla e, se necessario, corregge il separatore dopo il numero
    % (il controllo viene fatto sul nome originale, è solo informativo;
    %  la struttura ClCa usa fieldName come stringa, non modifica il file)
    afterField = regexp(fname, ['C', tok{1}, '([^0-9])'], 'tokens', 'once');
    if ~isempty(afterField)
        sep = afterField{1};
        if ~strcmp(sep, '-') && ~strcmp(sep, '_')
            % Il separatore non è né '-' né '_': situazione inattesa
            fprintf('  [INFO] Separatore inatteso dopo "%s" nel file: %s\n', ...
                fieldName, fname);
        end
    end
    % Nota: la correzione effettiva del nome file NON viene fatta qui;
    % il Field estratto è solo usato come metadato nella struttura ClCa.
end


% -------------------------------------------------------------------------
function newRow = buildClCaRow(sessionDate, sessionMouse, fieldName, ...
                               sessionFolder, fname830, fname920, s830, s920)
% BUILDCLCAROW  Costruisce una singola riga della struttura ClCa.
%
% Input:
%   sessionDate   - data della sessione (stringa, es. '20260420')
%   sessionMouse  - identificativo del mouse (es. 'Mouse17')
%   fieldName     - nome del campo (es. 'C1')
%   sessionFolder - nome della cartella della sessione
%   fname830      - nome del file a 830 nm
%   fname920      - nome del file a 920 nm
%   s830          - struttura out_struct del file a 830 nm
%   s920          - struttura out_struct del file a 920 nm
%
% Output:
%   newRow - struttura con tutti i campi di ClCa

    % Campi identificativi della riga (condivisi tra 830 e 920)
    newRow.Date     = sessionDate;
    newRow.Mouse    = sessionMouse;
    newRow.Field    = fieldName;
    newRow.Stimulus = '';          % da riempire in uno script successivo

    % Provenienza
    newRow.FoldName     = sessionFolder;
    newRow.FileName_830 = fname830;
    newRow.FileName_920 = fname920;

    % ------ Dati specifici per 830 nm ------------------------------------

    % G_sub_dark: canale verde con dark sottratto
    if isfield(s830, 'G_sub_dark')
        newRow.G_sub_dark_830 = s830.G_sub_dark;
    else
        newRow.G_sub_dark_830 = [];
        warning('Campo G_sub_dark assente in: %s', fname830);
    end

    % R_sub_dark: canale rosso con dark sottratto
    if isfield(s830, 'R_sub_dark')
        newRow.R_sub_dark_830 = s830.R_sub_dark;
    else
        newRow.R_sub_dark_830 = [];
        warning('Campo R_sub_dark assente in: %s', fname830);
    end

    % Time: vettore dei tempi
    if isfield(s830, 'Time')
        newRow.Time_830 = s830.Time;
    else
        newRow.Time_830 = [];
        warning('Campo Time assente in: %s', fname830);
    end

    % Stimulus_time: tempi degli stimoli
    if isfield(s830, 'Stimulus_time')
        newRow.Stimulus_time_830 = s830.Stimulus_time;
    else
        newRow.Stimulus_time_830 = [];
        warning('Campo Stimulus_time assente in: %s', fname830);
    end

    % ------ Dati specifici per 920 nm ------------------------------------

    if isfield(s920, 'G_sub_dark')
        newRow.G_sub_dark_920 = s920.G_sub_dark;
    else
        newRow.G_sub_dark_920 = [];
        warning('Campo G_sub_dark assente in: %s', fname920);
    end

    if isfield(s920, 'R_sub_dark')
        newRow.R_sub_dark_920 = s920.R_sub_dark;
    else
        newRow.R_sub_dark_920 = [];
        warning('Campo R_sub_dark assente in: %s', fname920);
    end

    if isfield(s920, 'Time')
        newRow.Time_920 = s920.Time;
    else
        newRow.Time_920 = [];
        warning('Campo Time assente in: %s', fname920);
    end

    if isfield(s920, 'Stimulus_time')
        newRow.Stimulus_time_920 = s920.Stimulus_time;
    else
        newRow.Stimulus_time_920 = [];
        warning('Campo Stimulus_time assente in: %s', fname920);
    end

end