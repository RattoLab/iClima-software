%%
%Second part of the analysis for the extraction of visual responses of
%different imaging trials in different animals (differet type of visual
%stimulations)

roilogi   = [final_struct.IsRoi]';
daylogi   = [final_struct.DayTime]'>=10 & [final_struct.DayTime]'<=14;
nightlogi = [final_struct.DayTime]'>=21 | [final_struct.DayTime]'<=3;
isoflogi  = [final_struct.IsIsoflu]';
mice      = unique({final_struct.Mouse});
isneurop  = strcmp({final_struct(:).Roi_type},'neurop')';
isIn      = strcmp({final_struct(:).Microscope},'CNR_IN')';
ibf       = arrayfun( @(x)string(x.Date),final_struct) > "20250200";
valori={final_struct.exclude};
valori(cellfun(@isempty,valori))={0};
[final_struct.exclude]=valori{:};
isempty(final_struct.exclude)=0;
exclude   = [final_struct.exclude]';

smooth = [3 0];

%% Select GRID0 long
stimlogi = cellfun(@(x)( ~isempty(regexp(x,'GRID\d+-10s-x1','once'))),...
                        {final_struct.StimName})';
t_sel = [-10 25];
t_bas = [-10  0];

%% Select GRID0 short
stimlogi = cellfun(@(x)( strcmp(x,'GRID0-5s-x1w')|...
                        strcmp(x,'GRID0-5s-x1b')|...
                        strcmp(x,'GRID0-5s-x1f')),{final_struct.StimName})';
t_sel = [-10 20];
t_bas = [-10  0];

%% STEP
stimlogi = cellfun(@(x)( ~isempty(regexpi(x,'FLASH','once'))),...
                        {final_struct.StimName})';
t_sel = [-10 35];
t_bas = [-10  0];

%% run 

sel = find(daylogi & isoflogi & ~isneurop & ibf & ~roilogi & stimlogi & ~exclude);
%sel = find(daylogi & ~isoflogi & ~isneurop & ibf & ~roilogi & stimlogi & ~exclude);
%sel = find(isoflogi & ~isneurop & ibf & ~roilogi & stimlogi & ~exclude);
%sel = find(~isoflogi & ~isneurop & ibf & ~roilogi & stimlogi & ~exclude);

names = "";

for i=1:numel(sel)
    t_temp  = final_struct(sel(i)).t_int;
    r  = movmean( 1./final_struct(sel(i)).ratio_int, smooth, 1);
    r0 = mean( r(t_temp>=t_bas(1) & t_temp<t_bas(2),:,:), [1 2], 'omitnan');
    x = (r-r0)./r0;
    xx = mean( x( t_temp>=t_sel(1) & t_temp<= t_sel(2),:,:), [2 3], 'omitnan');

    if i==1
        drr = nan(size(xx,1), numel(sel));
        drr(:,1) = xx;
        t = t_temp(t_temp>=t_sel(1) & t_temp<= t_sel(2));
    end
    
    drr(:,i) = xx*100;
    names(i) = string(final_struct(sel(i)).Mouse);
end


%% plot all sessions
figure
plot(t,drr,'color',[0.5 0.5 0.5])
hold on
plot(t,mean(drr,2,'omitnan'),'k','LineWidth',2)
ylim([-2 8])

%% plot all mice

micenames = unique(names);
figure
tiledlayout('flow')
avgs = nan(size(drr,1),numel(micenames));
for i=1:numel(micenames)
    now = find(strcmp(micenames(i),names));
    n = numel(now);
    tmp = drr(:,now);
    avgs(:,i) = mean(tmp,2,'omitnan');
    nexttile; hold on
    plot(t,tmp,'color',[0.5 0.5 0.5])
    plot(t,avgs(:,i),'k','LineWidth',2)
    title(sprintf('%s \nn sessions = %i',micenames(i),n))
    ylim([-2 8])
end
nexttile
plot(t,avgs)
hold on
plot(t,mean(avgs,2,'omitnan'),'k','LineWidth',3)
ylim([-2 8])

%%
risultato=0;
for i=1:numel(sel)

f=size(final_struct(sel(i)).G_sub_dark,3);
risultato=risultato+f;
end



%%
% 1. Filtro
target_mice = ["Mouse2", "Mouse6", "Mouse7", "Mouse8"];
mice_logi = ismember({final_struct.Mouse}, target_mice)';
sel_2 = find(daylogi & isoflogi & ~isneurop & ibf & ~roilogi & stimlogi & ~exclude & mice_logi);

% 2. Inizializzazione contenitore risultati
% Usiamo una mappa o una struct dinamica per salvare i dati per topo
for m = 1:numel(target_mice)
    results_by_mouse.(target_mice(m)).drr = [];
    results_by_mouse.(target_mice(m)).n_roi = 0;
    results_by_mouse.(target_mice(m)).n_sessions = 0;
end

% 3. Ciclo di analisi
for i = 1:numel(sel_2)
    idx = sel_2(i);
    mouse_name = string(final_struct(idx).Mouse);
    
    % Calcolo traccia (tua logica)
    t_temp  = final_struct(idx).t_int;
    r  = movmean(1./final_struct(idx).ratio_int, smooth, 1);
    r0 = mean(r(t_temp>=t_bas(1) & t_temp<t_bas(2),:,:), [1 2], 'omitnan');
    x = (r-r0)./r0;
    xx = mean(x(t_temp>=t_sel(1) & t_temp<= t_sel(2),:,:), [2 3], 'omitnan');
    
    % Aggiornamento contatore ROI
    dims = size(final_struct(idx).drr_int);
    n_roi_session = 1; if length(dims) >= 3, n_roi_session = dims(3); end
    
    % Salvataggio nei contenitori separati
    results_by_mouse.(mouse_name).drr = [results_by_mouse.(mouse_name).drr, xx*100];
    results_by_mouse.(mouse_name).n_roi = results_by_mouse.(mouse_name).n_roi + n_roi_session;
    results_by_mouse.(mouse_name).n_sessions = results_by_mouse.(mouse_name).n_sessions + 1;
end

% 4. Stampa dei risultati e Plot separato
figure; tiledlayout('flow');
for m = 1:numel(target_mice)
    name = target_mice(m);
    data = results_by_mouse.(name);
    
    fprintf('Topo %s: %d sessioni, %d ROI totali.\n', name, data.n_sessions, data.n_roi);
    
    % Plot individuale
    nexttile;
    plot(t, data.drr, 'color', [0.7 0.7 0.7]); hold on;
    plot(t, mean(data.drr, 2, 'omitnan'), 'r', 'LineWidth', 2);
    title(sprintf('%s (%d ROI)', name, data.n_roi));
    ylim([-2 8]); grid on;
end