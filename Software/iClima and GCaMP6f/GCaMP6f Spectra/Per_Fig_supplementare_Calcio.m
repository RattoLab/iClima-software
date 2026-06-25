%%For the spectra of the GCaMP6f 
% plot_ratios_wavelength.m
% Plotta automaticamente tutte le variabili "ratio*nm*" presenti nel
% workspace base contro il vettore tempo "t", ordinandole per
% lunghezza d'onda crescente. Genera sia una figura con subplot
% separati sia una figura con tutte le curve sovrapposte.
%
% Convenzioni assunte:
%   - le variabili si chiamano ratio<wavelength>nm[<suffisso opzionale>]
%     es: ratio800nm, ratio1000nm, ratio830nm1
%   - "t" è il vettore tempo (stessa lunghezza dei ratio*)
%   - variabili con lo stesso numero di wavelength ma suffisso diverso
%     (es. ratio830nm e ratio830nm1) vengono trattate come appartenenti
%     alla stessa wavelength e quindi plottate vicine/insieme

clear varNames wavelengths

%% 1. Trova il vettore tempo
if ~evalin('base', 'exist(''t'',''var'')')
    error('Variabile tempo "t" non trovata nel workspace base.');
end
t = evalin('base', 't');
t = t(:); % forza colonna

%% 2. Individua tutte le variabili che matchano il pattern ratio<numero>nm...
allVars = evalin('base', 'who');
pattern = '^ratio(\d+)nm(\d*)$';

varNames = {};
wavelengths = [];
suffixes = [];

for i = 1:numel(allVars)
    name = allVars{i};
    tok = regexp(name, pattern, 'tokens', 'once');
    if ~isempty(tok)
        varNames{end+1} = name; %#ok<AGROW>
        wavelengths(end+1) = str2double(tok{1}); %#ok<AGROW>
        if isempty(tok{2})
            suffixes(end+1) = 0; %#ok<AGROW>
        else
            suffixes(end+1) = str2double(tok{2}); %#ok<AGROW>
        end
    end
end

if isempty(varNames)
    error('Nessuna variabile del tipo ratio<wavelength>nm trovata nel workspace.');
end

%% 3. Ordina per wavelength crescente (e per suffisso a parità di wavelength)
[~, sortIdx] = sortrows([wavelengths(:), suffixes(:)], [1 2]);
varNames    = varNames(sortIdx);
wavelengths = wavelengths(sortIdx);
suffixes    = suffixes(sortIdx);

nVars = numel(varNames);

% Costruisce le label per le legende/titoli, distinguendo i duplicati
labels = cell(1, nVars);
for i = 1:nVars
    if suffixes(i) == 0
        labels{i} = sprintf('%d nm', wavelengths(i));
    else
        labels{i} = sprintf('%d nm (rep. %d)', wavelengths(i), suffixes(i));
    end
end

fprintf('Trovate %d variabili, ordine di plot:\n', nVars);
for i = 1:nVars
    fprintf('  %2d) %-15s -> %s\n', i, varNames{i}, labels{i});
end

%% 4. Carica i dati in una cella per riuso comodo
data = cell(1, nVars);
for i = 1:nVars
    v = evalin('base', varNames{i});
    data{i} = v(:); % forza colonna
end

%% 4bis. Smoothing causale (movmean "solo passato") con gestione NaN
% windowLen = numero di campioni della finestra causale (modificabile).
% Una finestra causale di N campioni usa, per il campione i, i valori
% da (i-N+1) a i, MAI valori futuri: questo è essenziale qui perché
% t=0 corrisponde allo stimolo, e non vogliamo che lo smoothing "veda"
% la risposta prima che lo stimolo sia avvenuto.
windowLen = 20; % <-- cambia qui per modificare la lunghezza della finestra

dataSm = cell(1, nVars);
for i = 1:nVars
    % movmean con finestra asimmetrica [windowLen-1, 0]:
    % usa solo il campione corrente e i (windowLen-1) precedenti.
    % 'omitnan' = i NaN in finestra vengono ignorati, la media è fatta
    % solo sui valori validi presenti (se la finestra è tutta NaN,
    % il risultato resta NaN).
    dataSm{i} = movmean(data{i}, [windowLen-1, 0], 'omitnan');

    % Salva la versione smoothed nel workspace base come nuova variabile
    smName = [varNames{i} '_sm'];
    assignin('base', smName, dataSm{i});
end

fprintf('\nCreate %d variabili smoothed (movmean causale, finestra = %d campioni):\n', nVars, windowLen);
for i = 1:nVars
    fprintf('  %s -> %s_sm\n', varNames{i}, varNames{i});
end

% Colormap coerente tra le due figure, costruita per massimizzare la
% distinguibilità tra curve adiacenti (importante con 10+ variabili,
% dove colormap "lisce" come parula/jet rendono simili i colori vicini
% in sequenza). Strategia: hue equidistanziato sul cerchio cromatico,
% con saturazione/luminosità alternate a step per aumentare ulteriormente
% il contrasto tra wavelength consecutive.
hues = (0:nVars-1)' / nVars;                 % hue equidistanziati [0,1)
sat  = repmat([0.85; 0.55], ceil(nVars/2), 1); sat = sat(1:nVars);   % alterna saturazione
val  = repmat([0.85; 0.65], ceil(nVars/2), 1); val = val(1:nVars);  % alterna luminosità
cmap = hsv2rgb([hues, sat, val]);

%% 5. FIGURA 1: subplot separati, uno per ciascuna wavelength (dati smoothed)
figure('Name', 'Ratio vs tempo - subplot (smoothed)', 'Color', 'w');
for i = 1:nVars
    subplot(nVars, 1, i);
    plot(t, dataSm{i}, 'Color', cmap(i,:), 'LineWidth', 1.2);
    ylabel(labels{i}, 'FontSize', 8);
    grid on;
    if i < nVars
        set(gca, 'XTickLabel', []);
    else
        xlabel('Tempo');
    end
end
sgtitle(sprintf('Ratio per ciascuna lunghezza d''onda (ordine crescente) - movmean causale, finestra %d campioni', windowLen));

%% 6. FIGURA 2: tutte le curve sovrapposte nello stesso plot (dati smoothed)
figure('Name', 'Ratio vs tempo - overlay (smoothed)', 'Color', 'w');
hold on;
for i = 1:nVars
    plot(t, dataSm{i}, 'Color', cmap(i,:), 'LineWidth', 1.2, ...
        'DisplayName', labels{i});
end
hold off;
grid on;
xlabel('Tempo');
ylabel('Ratio (smoothed)');
title(sprintf('Ratio vs tempo, tutte le lunghezze d''onda (crescente) - movmean causale, finestra %d campioni', windowLen));
legend('show', 'Location', 'bestoutside');