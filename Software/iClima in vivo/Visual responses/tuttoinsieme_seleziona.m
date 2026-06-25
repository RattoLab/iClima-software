%%
%First part of the analysis for the extraction of visual responses of
%different imaging trials in different animals (differet type of visual
%stimulations)
tic
disp('Loading...')
% load('\\146.48.73.104\Data\Data_QNAP\Chloride Imaging\AAV iClima\IBF\Analyzed\MERGED\19-Aug-2025 17-00-44Merged_structure.mat');
% load('\\NAS617510\Data\Data_QNAP\Chloride Imaging\AAV iClima\IBF\Analyzed\MERGED\all_mice_20250915.mat');
load("\\146.48.73.104\Data\Data_QNAP\Paper Cl1\Visual Responses for paper\visual Responses_iClima\Dataset\all_mice_20250915.mat");
% the loaded structure is called "final_struct".
fprintf('loaded in %.2f s.\n', toc)
filename = {final_struct.Filename};
date = {final_struct.Date};
dayt = daytime_str2double( {final_struct.Day_time} );

stimnames = load("stimnames.mat");
[stn, sti] = assign_names(filename,date,stimnames);

% we have some duplicated data
duplic_idx = [...
    36 44 47 51 53 54 59 70 104 108 112 114 115 120 123 127 129 130 135 136 137 139 146 173 180 227 228 ... % flash
    55 56 57 58 60 116 117 118 119 121 131 132 133 134 136 138 140 141 143 144 147 148 150 151 153 154 156 157 158 161 162 164 165 166 167 168 171 195 365 489 523 529 574 581 ... % grid
    142 145 149 152 155 159 160 172 174 181 581 ... % grid 10s
    ];

% to excl by giacomo
excl_giack = [365 489 523 529 574 581];

for i=1:numel(final_struct)
   
    % assign official names
    final_struct(i).StimName = stn{i}; %#ok<*SAGROW>
    final_struct(i).Stim_i = sti(i);

    % add the "IsRoi" field
    final_struct(i).IsRoi = ~isempty( regexpi( filename{i}, 'ROI', 'once')) |...
            1/mean(diff(final_struct(i).Time)) > 20; % faster than 20 Hz

    % add the "Time" (numeric) field
    final_struct(i).DayTime = dayt(i);

    % add the "IsIsoflu" field
    final_struct(i).IsIsoflu = ~isempty( regexpi( filename{i}, 'isoon', 'once')) |...
            isequal(final_struct(i).isoflurane,1);

    % add the "IsNeuropile" field
    final_struct(i).IsNeuropile = ~isempty( regexpi( filename{i}, 'neuropile', 'once'));

    % Clean the "Field" field. It should be identified by the letter C
    % (campo) followed by a progressive number and optionally by the string
    % "ROI". If the character "-" is found, discard anything that comes 
    % after it, dash included.
    dash = regexp(final_struct(i).Field,"-");
    if ~isempty(dash)
        final_struct(i).Field = final_struct(i).Field(1:dash(1)-1);
    end

    % add the field exclude
    final_struct(i).exclude = false;

    % exclude the duplicated data
    if ismember(i,duplic_idx)
        final_struct(i).exclude = true;
    end
    
    % excluded by giak
    if ismember(i,excl_giack)
        final_struct(i).exclude = true;
    end
end



%% recompute

disp('Recomputing deltaR/R...')

% Select a time time window fot the computation of R0 (baseline)
rat0_window = [-10 0];

% Do you want to compute the R0 on the average of all trials?
rat0_on_avg = false;
% rat0_on_avg = true;


% If you compute the R0 on the average of all the trials, exclude from the
% average the those trials when the mouse was running
exclude_movement_for_r0 = false; % it applies only if rat0_on_avg is true!
% exclude_movement_for_r0 = true; % it applies only if rat0_on_avg is true!

speed_thr = 0.02; % Threshold for run (m/s)

% Apply BT correction (BT values at 830nm)
alpha = 0.241; % use ZERO if you want no BT subtraction
beta = 0.0839;
% alpha=0;
% beta=0;

for i=1:numel(final_struct)

    s = final_struct(i);
    t = s.Time;
    % remove bt
    g = s.G_sub_dark - s.R_sub_dark.*alpha;
    r = s.R_sub_dark - s.G_sub_dark.*beta;
    
     ratio = g./r; % green / red ratio
    %ratio = r./g; % red / green ratio
    try
        if rat0_on_avg
            if exclude_movement_for_r0
                run = s.Speed;
                running = abs( run( t>-rat0_window(1) & t<rat0_window(2), :) )...
                    < speed_thr;
                trials_for_r0 = all( running, 1);
            else
                trials_for_r0 = true( 1, size( ratio, 2));
            end

            if trials_for_r0>0
                % average trials
                avg_trial = mean( ratio(:,trials_for_r0,:), 2, 'omitnan');
                % compute baseline ratio. This is R0
                rat0 = mean( avg_trial( t<rat0_window(2) & t>rat0_window(1), :, :), 1);
                % compute delta R over R0. Same as calcium imaging
                drat_over_rat = ( ratio - rat0) ./ rat0;
                drat_over_rat_trials = squeeze( mean( drat_over_rat, 2, 'omitnan'));
                avg_response = mean( drat_over_rat_trials, 2, 'omitnan');
            else
                avg_trial = [];
                rat0 = [];
                drat_over_rat = [];
                drat_over_rat_trials=[];
                avg_response=[];
            end
        else
            % compute baseline ratio. This is R0
            rat0 = mean( ratio( t < 0 & t>-10, :, :), 1);
            % compute delta R over R0. Same as calcium imaging
            drat_over_rat = ( ratio - rat0) ./ rat0;
            % average over trials
            drat_over_rat_trials = squeeze( mean( drat_over_rat, 2, 'omitnan'));
            avg_response = mean( drat_over_rat_trials, 2, 'omitnan');
        end
        final_struct(i).r0_recomp = rat0;
        final_struct(i).drr       = drat_over_rat;
        final_struct(i).drr_avg   = drat_over_rat_trials;
        final_struct(i).Ratio     = ratio; 
    catch ERR
        final_struct(i).exclude = false;
        warning('\n\nError while recomputing trace %i',i)
        warning(getReport(ERR));
    end
end
disp("DeltaR/R recomputed.")



%% align and resample
time_from = -30;
time_to   = 90;
disp("Aligning and resampling...")
for i=1:numel(final_struct)
    try
        current_sp = median(diff(final_struct(i).Time));
        if current_sp>0.2 % <5Hz (normal frame rate 2 Hz)
            % current_sp = floor(current_sp*10)/10; % round to 100ms precision
            current_sp = 0.5; % 2 Hz
        else % >5Hz (frame rate  for rois = 30 Hz = 0.033 s)
            % current_sp = floor(current_sp*1000)/1000; % truncates to 10ms precision
            current_sp = 0.033; % circa 30 Hz
        end
        before_n = round(time_from/current_sp);
        after_n  = round(time_to  /current_sp);
        newt = (before_n*current_sp:current_sp:after_n*current_sp)'; % 2 Hz
        nt = numel(newt);

        current_ratio = final_struct(i).Ratio;
        current_drr = final_struct(i).drr;
        [ ~, ntrials, nrois] = size(current_drr);

        tmprat= zeros(nt, ntrials, nrois);
        tmpr  = zeros(nt, ntrials, nrois);
        tmpv  = zeros(nt,ntrials);
        for itr=1:ntrials
            for iroi=1:nrois
                tmprat(:,itr,iroi) = interp1( final_struct(i).Time, ...
                    current_ratio(:, itr, iroi), ...
                    newt,...
                    'linear',nan);
                tmpr(:, itr,iroi) = interp1( final_struct(i).Time, ...
                    current_drr(:, itr, iroi), ...
                    newt,...
                    'linear',nan);
            end
            tmpv(:, itr) = interp1( final_struct(i).Time, ...
                final_struct(i).Speed(:,itr), ...
                newt,...
                'linear',nan);
        end
        final_struct(i).t_int       = newt;
        final_struct(i).ratio_int   = tmprat;
        final_struct(i).drr_int     = tmpr;
        final_struct(i).drr_avg_int = squeeze( mean( final_struct(i).drr_int, 2, 'omitnan'));
        final_struct(i).speed_int   = tmpv;
    catch ERR
        final_struct(i).exclude = true;
        warning('\n\nError while resampling trace %i',i)
        warning(getReport(ERR));
    end
end
disp("Aligned and resampled.")


%% labels

roilogi   = [final_struct.IsRoi]';
daylogi   = [final_struct.DayTime]'>=10 & [final_struct.DayTime]'<=14;
nightlogi = [final_struct.DayTime]'>=21 | [final_struct.DayTime]'<=3;
isoflogi  = [final_struct.IsIsoflu]';
mice      = unique({final_struct.Mouse});
isneurop  = [final_struct.IsNeuropile]';
isIn      = strcmp({final_struct(:).Microscope},'CNR_IN')';

%% Select GRID0 SHORT

grid0logi = cellfun(@(x)( strcmp(x,'GRID0-5s-x1w')|...
                        strcmp(x,'GRID0-5s-x1b')|...
                        strcmp(x,'GRID0-5s-x1f')),{final_struct.StimName})';

% display how many
for i=1:numel(mice)
    thismousei = strcmp({final_struct.Mouse},mice{i})';
    fprintf('%s:\n',mice{i})
    fprintf('Day: %i awake, %i isoflu\n',...
        sum( thismousei & daylogi & ~isoflogi & grid0logi & ~roilogi),...
        sum( thismousei & daylogi &  isoflogi & grid0logi & ~roilogi))
    fprintf('Night: %i awake, %i isoflu\n\n',...
        sum( thismousei & nightlogi & ~isoflogi & grid0logi & ~roilogi),...
        sum( thismousei & nightlogi &  isoflogi & grid0logi & ~roilogi))
end

% Plot 
for i=1:numel(mice)
    thismousei = strcmp({final_struct.Mouse},mice{i})';
    % fprintf('%s:\n',mice{i})
    dayi      = find(thismousei & grid0logi & not(roilogi) & daylogi   & ~isoflogi);
    dayisoi   = find(thismousei & grid0logi & not(roilogi) & daylogi   &  isoflogi);
    nighti    = find(thismousei & grid0logi & not(roilogi) & nightlogi & ~isoflogi);
    nightisoi = find(thismousei & grid0logi & not(roilogi) & nightlogi &  isoflogi);
    % fprintf('%s:\n',mice{i})
    % disp(cat(1,dayi,dayisoi,nighti,nightisoi))
    figure('Name',mice{i})
    tiledlayout(1,2)
    nexttile
    hold on
    plot([-40 0 0 5 5 100],[0 0 1 1 0 0],'r')
    for ii=1:numel(dayi)
        plot( final_struct(dayi(ii)).t_int,...
            mean( final_struct(dayi(ii)).drr_avg_int.*100, 2, 'omitnan'),'g');
    end
    for ii=1:numel(nighti)
        plot( final_struct(nighti(ii)).t_int,...
            mean( final_struct(nighti(ii)).drr_avg_int.*100, 2, 'omitnan'),'b');
    end
    ylabel('\Deltaratio/ratio_0 (%)')
    xlabel('Time (s)')
    set(gca,'ydir','reverse','ylim',[-10 8],'xlim',[-40 100])
    title('Awake')

    nexttile
    hold on
    plot([-40 0 0 5 5 100],[0 0 1 1 0 0],'r')
    for ii=1:numel(dayisoi)
        plot( final_struct(dayisoi(ii)).t_int,...
            mean( final_struct(dayisoi(ii)).drr_avg_int.*100, 2, 'omitnan'),'g');
    end
    for ii=1:numel(nightisoi)
        plot( final_struct(nightisoi(ii)).t_int,...
            mean( final_struct(nightisoi(ii)).drr_avg_int.*100, 2, 'omitnan'),'b');
    end
    ylabel('\Deltaratio/ratio_0 (%)')
    xlabel('Time (s)')
    set(gca,'ydir','reverse','ylim',[-10 8],'xlim',[-40 100])
    title('Isoflurane')
end


%% Select GRID "DEVIANT"

stimlogi = cellfun(@(x)( ~isempty(regexp(x,'GRID\d+-\d+gf-5s-x1','once'))),...
                        {final_struct.StimName})';

% display how many
for i=1:numel(mice)
    thismousei = strcmp({final_struct.Mouse},mice{i})';
    fprintf('%s:\n',mice{i})
    fprintf('Day: %i awake, %i isoflu\n',...
        sum( thismousei & daylogi & ~isoflogi & stimlogi & ~roilogi),...
        sum( thismousei & daylogi &  isoflogi & stimlogi & ~roilogi))
    fprintf('Night: %i awake, %i isoflu\n\n',...
        sum( thismousei & nightlogi & ~isoflogi & stimlogi & ~roilogi),...
        sum( thismousei & nightlogi &  isoflogi & stimlogi & ~roilogi))
end

% Plot 

for i=1:numel(mice)
    thismousei = strcmp({final_struct.Mouse},mice{i})';
    % fprintf('%s:\n',mice{i})
    dayi      = find(thismousei & stimlogi & not(roilogi) & daylogi   & ~isoflogi);
    dayisoi   = find(thismousei & stimlogi & not(roilogi) & daylogi   &  isoflogi);
    nighti    = find(thismousei & stimlogi & not(roilogi) & nightlogi & ~isoflogi);
    nightisoi = find(thismousei & stimlogi & not(roilogi) & nightlogi &  isoflogi);
    % fprintf('%s:\n',mice{i})
    % disp(cat(1,dayi,dayisoi,nighti,nightisoi))
    figure('Name',mice{i})
    tiledlayout(1,2)
    nexttile
    hold on
    plot([-40 0 0 5 5 10 10 100],[0 0 1 1 -1 -1 0 0],'r')
    for ii=1:numel(dayi)
        plot( final_struct(dayi(ii)).t_int,...
            mean( final_struct(dayi(ii)).drr_avg_int.*100, 2, 'omitnan'),'g');
    end
    for ii=1:numel(nighti)
        plot( final_struct(nighti(ii)).t_int,...
            mean( final_struct(nighti(ii)).drr_avg_int.*100, 2, 'omitnan'),'b');
    end
    ylabel('\Deltaratio/ratio_0 (%)')
    xlabel('Time (s)')
    set(gca,'ydir','reverse','ylim',[-10 8],'xlim',[-40 100])
    title('Awake')

    nexttile
    hold on
    plot([-40 0 0 5 5 10 10 100],[0 0 1 1 -1 -1 0 0],'r')
    for ii=1:numel(dayisoi)
        plot( final_struct(dayisoi(ii)).t_int,...
            mean( final_struct(dayisoi(ii)).drr_avg_int.*100, 2, 'omitnan'),'g');
    end
    for ii=1:numel(nightisoi)
        plot( final_struct(nightisoi(ii)).t_int,...
            mean( final_struct(nightisoi(ii)).drr_avg_int.*100, 2, 'omitnan'),'b');
    end
    ylabel('\Deltaratio/ratio_0 (%)')
    xlabel('Time (s)')
    set(gca,'ydir','reverse','ylim',[-10 8],'xlim',[-40 100])
    title('Isoflurane')
end

%% Select GRID LONG

stimlogi = cellfun(@(x)( ~isempty(regexp(x,'GRID\d+-10s-x1','once'))),...
                        {final_struct.StimName})';

% display how many
for i=1:numel(mice)
    thismousei = strcmp({final_struct.Mouse},mice{i})';
    fprintf('%s:\n',mice{i})
    fprintf('Day: %i awake, %i isoflu\n',...
        sum( thismousei & daylogi & ~isoflogi & stimlogi & ~roilogi),...
        sum( thismousei & daylogi &  isoflogi & stimlogi & ~roilogi))
    fprintf('Night: %i awake, %i isoflu\n\n',...
        sum( thismousei & nightlogi & ~isoflogi & stimlogi & ~roilogi),...
        sum( thismousei & nightlogi &  isoflogi & stimlogi & ~roilogi))
end

% Plot 

for i=1:numel(mice)
    thismousei = strcmp({final_struct.Mouse},mice{i})';
    % fprintf('%s:\n',mice{i})
    dayi      = find(thismousei & stimlogi & not(roilogi) & daylogi   & ~isoflogi);
    dayisoi   = find(thismousei & stimlogi & not(roilogi) & daylogi   &  isoflogi);
    nighti    = find(thismousei & stimlogi & not(roilogi) & nightlogi & ~isoflogi);
    nightisoi = find(thismousei & stimlogi & not(roilogi) & nightlogi &  isoflogi);
    % fprintf('%s:\n',mice{i})
    % disp(cat(1,dayi,dayisoi,nighti,nightisoi))
    figure('Name',mice{i})
    tiledlayout(1,2)
    nexttile
    hold on
    plot([-40 0 0 10 10 100],[0 0 1 1 0 0],'r')
    for ii=1:numel(dayi)
        plot( final_struct(dayi(ii)).t_int,...
            mean( final_struct(dayi(ii)).drr_avg_int.*100, 2, 'omitnan'),'g');
    end
    for ii=1:numel(nighti)
        plot( final_struct(nighti(ii)).t_int,...
            mean( final_struct(nighti(ii)).drr_avg_int.*100, 2, 'omitnan'),'b');
    end
    ylabel('\Deltaratio/ratio_0 (%)')
    xlabel('Time (s)')
    set(gca,'ydir','reverse','ylim',[-10 8],'xlim',[-40 100])
    title('Awake')

    nexttile
    hold on
    plot([-40 0 0 10 10 100],[0 0 1 1 0 0],'r')
    for ii=1:numel(dayisoi)
        plot( final_struct(dayisoi(ii)).t_int,...
            mean( final_struct(dayisoi(ii)).drr_avg_int.*100, 2, 'omitnan'),'g');
    end
    for ii=1:numel(nightisoi)
        plot( final_struct(nightisoi(ii)).t_int,...
            mean( final_struct(nightisoi(ii)).drr_avg_int.*100, 2, 'omitnan'),'b');
    end
    set(gca,'ydir','reverse','ylim',[-10 8],'xlim',[-40 100])
    ylabel('\Deltaratio/ratio_0 (%)')
    xlabel('Time (s)')
    title('Isoflurane')
end

%% Select FLASH

stimlogi = cellfun(@(x)( ~isempty(regexpi(x,'FLASH','once'))),...
                        {final_struct.StimName})';

% display how many
for i=1:numel(mice)
    thismousei = strcmp({final_struct.Mouse},mice{i})';
    fprintf('%s:\n',mice{i})
    fprintf('Day: %i awake, %i isoflu\n',...
        sum( thismousei & daylogi & ~isoflogi & stimlogi & ~roilogi),...
        sum( thismousei & daylogi &  isoflogi & stimlogi & ~roilogi))
    fprintf('Night: %i awake, %i isoflu\n\n',...
        sum( thismousei & nightlogi & ~isoflogi & stimlogi & ~roilogi),...
        sum( thismousei & nightlogi &  isoflogi & stimlogi & ~roilogi))
end

% Plot 

for i=1:numel(mice)
    thismousei = strcmp({final_struct.Mouse},mice{i})';
    % fprintf('%s:\n',mice{i})
    dayi      = find(thismousei & stimlogi & not(roilogi) & daylogi   & ~isoflogi);
    dayisoi   = find(thismousei & stimlogi & not(roilogi) & daylogi   &  isoflogi);
    nighti    = find(thismousei & stimlogi & not(roilogi) & nightlogi & ~isoflogi);
    nightisoi = find(thismousei & stimlogi & not(roilogi) & nightlogi &  isoflogi);
    % fprintf('%s:\n',mice{i})
    % disp(cat(1,dayi,dayisoi,nighti,nightisoi))
    figure('Name',mice{i})
    tiledlayout(1,2)
    nexttile
    hold on
    plot([-40 0 0 20 20 100],[0 0 1 1 0 0],'r')
    for ii=1:numel(dayi)
        plot( final_struct(dayi(ii)).t_int,...
            mean( final_struct(dayi(ii)).drr_avg_int.*100, 2, 'omitnan'),'g');
    end
    for ii=1:numel(nighti)
        plot( final_struct(nighti(ii)).t_int,...
            mean( final_struct(nighti(ii)).drr_avg_int.*100, 2, 'omitnan'),'b');
    end
    ylabel('\Deltaratio/ratio_0 (%)')
    xlabel('Time (s)')
    set(gca,'ydir','reverse','ylim',[-10 8],'xlim',[-40 100])
    title('Awake')

    nexttile
    hold on
    plot([-40 0 0 20 20 100],[0 0 1 1 0 0],'r')
    for ii=1:numel(dayisoi)
        plot( final_struct(dayisoi(ii)).t_int,...
            mean( final_struct(dayisoi(ii)).drr_avg_int.*100, 2, 'omitnan'),'g');
    end
    for ii=1:numel(nightisoi)
        plot( final_struct(nightisoi(ii)).t_int,...
            mean( final_struct(nightisoi(ii)).drr_avg_int.*100, 2, 'omitnan'),'b');
    end
    set(gca,'ydir','reverse','ylim',[-10 8],'xlim',[-40 100])
    ylabel('\Deltaratio/ratio_0 (%)')
    xlabel('Time (s)')
    title('Isoflurane')
end


