%% Script fot the comparison of iClima and ClopHensor responses
%%select a directory with the analyzed mat files for either iClima or ClopHensor
clc; clear;
dataDir = uigetdir(pwd, 'Select folder containing your .mat files');
cd(dataDir);
files = dir('*.mat');
%% 

 %alpha= 0.241;  %sono gli alpha e beta vecchi, ovvero i più conservativi
 %beta= 0.0839 ;
 alpha= 0;  %sono gli alpha e beta vecchi, ovvero i più conservativi
 beta= 0 ;

summary = struct();
timeMask_boundaries = [-15 15];  % common time window
timeVecCommon = [];
% --- Load & preprocess data 
for i = 1:length(files)
    fname = files(i).name;
    disp(['Processing ', fname]);
    data = load(fname);
    s = struct2cell(data);
    s = s{1};
    summary(i).Mouse       = s.Mouse;
    summary(i).Date        = s.Date;
    summary(i).Time        = s.Time;
    summary(i).NumROIs     = length(s.Roi_labels);
    summary(i).DateTime    = datetime(s.Date, 'InputFormat', 'yyyyMMdd');
     % dato che erano stati acquisiti al 2P di IN canale verde e canale
    % rosso sono invertiti
    summary(i).G_sub_dark  = s.R_sub_dark;
    summary(i).R_sub_dark  = s.G_sub_dark;
    %summary(i).Ratio       = s.Ratio;
    g_BT_correct= (s.R_sub_dark-(alpha*s.G_sub_dark))/(1+alpha*beta);
    r_BT_correct= (s.G_sub_dark-(beta*s.R_sub_dark))/(1+alpha*beta); 
    summary(i).g_BT_correct= g_BT_correct;
    summary(i).r_BT_correct= r_BT_correct;
    %Ratio_true=g_BT_correct./r_BT_correct;
    Ratio_true=r_BT_correct./g_BT_correct;  %Rosso su verde
    summary(i).Ratio_true= Ratio_true;
   

   
    time = s.Time;
    timeMask = time >= timeMask_boundaries(1) & time < timeMask_boundaries(2);
    summary(i).timeMask = timeMask;

    % Find common time vector intersection
    if isempty(timeVecCommon)
        timeVecCommon = time(timeMask);
    else
        timeVecCommon = intersect(timeVecCommon, time(timeMask));
    end

    ratio = Ratio_true;  % [time x trials x ROIs]
      baselineMask = time >= -10 & time <= 0;
    rat0 = mean(ratio(baselineMask, :, :), 1);
    drat_over_rat = (ratio - rat0) ./ rat0;
    ratio_avg = movmean(100 * drat_over_rat,[3,0]); %for smooth filtering
    summary(i).ratio_avg = ratio_avg;  %  smoothed, all trials included
end




%%  Individual trials with average across rois
% --- Crop to common time points and concatenate trials
all_trials_concat = [];  % [time x all_trials]

for i = 1:length(summary)
    time_i = summary(i).Time;                      % time vector of this file
    ratio_avg_i_full = summary(i).ratio_avg;       % [time x trials x ROIs]

    nTrials = size(ratio_avg_i_full, 2);
    tmp_mat = NaN(length(timeVecCommon), nTrials); % [common_time x trials]

    for t = 1:nTrials
        % Extract trial data: [time x ROIs]
        trial_data = squeeze(ratio_avg_i_full(:, t, :));  % [time x ROIs]
        
        % Average over ROIs: [time x 1]
        trial_mean = mean(trial_data, 2, 'omitnan');

        % Interpolate to common time vector
        tmp_mat(:, t) = interp1(time_i, trial_mean, timeVecCommon, 'linear', NaN);
    end
    % Concatenate to full matrix
    all_trials_concat = [all_trials_concat, tmp_mat];  % [common_time x all_trials]
end

% --- Compute mean and SEM over all trials (no filtering)
mean_trace = mean(all_trials_concat, 2, 'omitnan');
sem_trace = std(all_trials_concat, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(all_trials_concat(1,:))));

% Ensure column vectors
timeVecCommon = timeVecCommon(:);
mean_trace = mean_trace(:);
sem_trace = sem_trace(:);

figure('Color', 'w'); hold on;
fill([timeVecCommon; flipud(timeVecCommon)], ...
     [mean_trace - sem_trace; flipud(mean_trace + sem_trace)], ...
     [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(timeVecCommon, mean_trace, 'b-', 'LineWidth', 2);
ylim([-2 8])
xlabel('Time (s)');
ylabel('ΔR/R (%)');
title('Average Ratio Response');
% grid on;
hold off;

%% Plot individual trials along with average trace and SEM

figure('Color', 'w'); hold on;

plot(timeVecCommon, all_trials_concat, 'Color', [0.7 0.7 0.7 0.3]);
fill([timeVecCommon; flipud(timeVecCommon)], ...
     [mean_trace - sem_trace; flipud(mean_trace + sem_trace)], ...
     [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(timeVecCommon, mean_trace, 'b-', 'LineWidth', 2);

num_trials = size(all_trials_concat, 2);

ylim([-2 8])
xlabel('Time (s)');
ylabel('ΔR/R (%)');
title(sprintf('Individual Trials (n = %d) with Average ', num_trials));
hold off;

%% Average each file across ROIs and trials, plot all file averages with overall average with SEM

file_avg_traces = NaN(length(timeVecCommon), length(summary));
figure('Color', 'w'); hold on;

for i = 1:length(summary)
    time_i = summary(i).Time;
    ratio_avg_i_full = summary(i).ratio_avg;  % [time x trials x ROIs]
    % Average over ROIs and trials for this file
    avg_trace_full = mean(mean(ratio_avg_i_full, 3, 'omitnan'), 2, 'omitnan');  % [time x 1]
    % Interpolate onto common time vector using linear method
    avg_trace_common = interp1(time_i, avg_trace_full, timeVecCommon, 'linear', NaN);
    % Store for overall average
    file_avg_traces(:, i) = avg_trace_common;
    % Plot this file in light grey
    plot(timeVecCommon, avg_trace_common, 'Color', [0.7 0.7 0.7 0.6]);
    picco(i)=max(avg_trace_common(timeVecCommon>0 & timeVecCommon<=10));  %serve per calcolar
end

% Compute overall mean and SEM across files (omit NaNs)
overall_mean = mean(file_avg_traces, 2, 'omitnan');
overall_sem = std(file_avg_traces, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(file_avg_traces), 2));

% Plot SEM shaded area
fill([timeVecCommon; flipud(timeVecCommon)], ...
     [overall_mean - overall_sem; flipud(overall_mean + overall_sem)], ...
     [0.9 0.6 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

plot(timeVecCommon, overall_mean, 'k-', 'LineWidth', 3);

% errorbar(timeVecCommon, overall_mean, overall_sem, 'LineStyle', 'none', 'Color', [0.4 0.6 0.1]);
ylim([-2 8])
xlabel('Time (s)');
ylabel('ΔR/R (%)');
nFiles = length(summary);
title(sprintf('Session (n = %d) averages (grey) and their mean (black) with SEM', nFiles));
% grid on;
hold off;

%% Plot file averages and mouse average + SEM for each mouse in the same figure

all_mice = {summary.Mouse};
unique_mice = unique(all_mice);
nMice = length(unique_mice);
colors = lines(nMice); % Distinct color for each mouse

figure('Color', 'w'); hold on;
legendHandles = gobjects(nMice + 1, 1); % +1 for grand mean
legendEntries = cell(nMice + 1, 1);     % corresponding labels
% Store mouse means to compute grand average later
mouse_means_all = NaN(length(timeVecCommon), nMice);
for m = 1:nMice
    mouse_id = unique_mice{m};
    idx_mouse = find(strcmp(all_mice, mouse_id));
    % Preallocate for session traces of this mouse
    file_avg_traces_mouse = NaN(length(timeVecCommon), length(idx_mouse));
    
    for j = 1:length(idx_mouse)
        i = idx_mouse(j);
        time_i = summary(i).Time;
        ratio_avg_i_full = summary(i).ratio_avg;  % [time x trials x ROIs]

        % Average over ROIs and trials
        avg_trace_full = mean(mean(ratio_avg_i_full, 3, 'omitnan'), 2, 'omitnan');

        % Interpolate to common time base
        avg_trace_common = interp1(time_i, avg_trace_full, timeVecCommon, 'linear', NaN);

        file_avg_traces_mouse(:, j) = avg_trace_common;

% %         OPTIONAL: Plot individual sessions (commented out for clarity)
%         plot(timeVecCommon, avg_trace_common, 'Color', [colors(m,:), 0.3]);
    end

    % Compute mean and SEM for this mouse
    mouse_mean = mean(file_avg_traces_mouse, 2, 'omitnan');
    mouse_means_all(:, m) = mouse_mean;  % store for grand average

    % Plot thick mouse mean
    h = plot(timeVecCommon, mouse_mean, '-', 'Color', colors(m,:), 'LineWidth', 1);
    legendHandles(m) = h;
    legendEntries{m} = sprintf('%s (n = %d)', mouse_id, length(idx_mouse));
end

% ---- Plot grand average (mean of mouse means) ----
grand_mean = mean(mouse_means_all, 2, 'omitnan');
h_grand = plot(timeVecCommon, grand_mean, 'm-', 'LineWidth', 3); % black dashed line
legendHandles(end) = h_grand;
legendEntries{end} = 'Mice Mean';

% ---- Final plot settings ----
ylim([-2 8]);
xlabel('Time (s)');
ylabel('ΔR/R (%)');
title('Averages of session averages for each mouse');

legend(legendHandles, legendEntries, 'Location', 'best');
hold off;

%% ---- Separate plots for each mouse ----

figure('Color', 'w');
nCols = ceil(sqrt(nMice)); % Square-ish layout
nRows = ceil(nMice / nCols);

for m = 1:nMice
    mouse_id = unique_mice{m};
    idx_mouse = find(strcmp(all_mice, mouse_id));
    
    file_avg_traces_mouse = NaN(length(timeVecCommon), length(idx_mouse));

    subplot(nRows, nCols, m); hold on;

    for j = 1:length(idx_mouse)
        i = idx_mouse(j);
        time_i = summary(i).Time;
        ratio_avg_i_full = summary(i).ratio_avg;

        avg_trace_full = mean(mean(ratio_avg_i_full, 3, 'omitnan'), 2, 'omitnan');
        avg_trace_common = interp1(time_i, avg_trace_full, timeVecCommon, 'linear', NaN);

        file_avg_traces_mouse(:, j) = avg_trace_common;

        % Plot individual file trace
        plot(timeVecCommon, avg_trace_common, 'Color', [colors(m,:), 0.3]);
    end

    % Plot mouse average
    mouse_mean = mean(file_avg_traces_mouse, 2, 'omitnan');
    plot(timeVecCommon, mouse_mean, '-', 'Color', colors(m,:), 'LineWidth', 2);
    
    ylim([-2 8]);
    xlim([min(timeVecCommon), max(timeVecCommon)]);
    
    % Add mouse ID and number of files
    title(sprintf('%s (n = %d)', mouse_id, length(idx_mouse)));

    if m > (nRows - 1) * nCols
        xlabel('Time (s)');
    end
    if mod(m - 1, nCols) == 0
        ylabel('ΔR/R (%)');
    end
    box off;
end

sgtitle(sprintf('Per-Mouse Session Averages\n(n is the session number)'));

%% Plot Grand Mean, Overall Mean, and Trial Mean Together
figure('Color', 'w'); hold on;

% Plot Trial Mean (blue solid)
plot(timeVecCommon, mean_trace, 'b-', 'LineWidth', 2);

% Plot Session Mean (black solid)
plot(timeVecCommon, overall_mean, 'k-', 'LineWidth', 2);

% Plot Grand Mean (mouse-level) (magenta dashed)
plot(timeVecCommon, grand_mean, 'm-', 'LineWidth', 2);

ylim([-1 4]);
xlim([-10 15]);
xlabel('Time (s)');
ylabel('ΔR/R (%)');
title('Comparison of Mean Traces Across Levels');

legend({'Trials Mean', 'Sessions Mean', 'Mice Mean'}, 'Location', 'best');
hold off;

%% Plot sessions per mouse sorted by DATE

for m = 1:nMice
    mouse_id = unique_mice{m};
    idx_mouse = find(strcmp(all_mice, mouse_id));
    
    % Sort by Date (string), not full datetime
    session_dates = {summary(idx_mouse).Date};  % Cell array of date strings
    [~, sort_idx] = sort(session_dates);        % Sort by string date
    sorted_idxs = idx_mouse(sort_idx);
    nSessions = length(sorted_idxs);

    % Set up figure
    figure('Name', ['Sessions for Mouse ' mouse_id], 'Color', 'w');
    nCols = ceil(sqrt(nSessions));
    nRows = ceil(nSessions / nCols);

    for j = 1:nSessions
        i = sorted_idxs(j);
        time_i = summary(i).Time;
        ratio_avg_i_full = summary(i).ratio_avg;

        % Average over ROIs and trials
        avg_trace_full = mean(mean(ratio_avg_i_full, 3, 'omitnan'), 2, 'omitnan');
        avg_trace_common = interp1(time_i, avg_trace_full, timeVecCommon, 'linear', NaN);

        subplot(nRows, nCols, j); hold on;
        plot(timeVecCommon, avg_trace_common, '-', 'Color', colors(m,:), 'LineWidth', 2);
        ylim([-2 8]);
        xlim([min(timeVecCommon), max(timeVecCommon)]);

        % Use Date as title (not full datetime)
        title(summary(i).Date, 'Interpreter', 'none', 'FontSize', 8);

        if j > (nRows - 1) * nCols
            xlabel('Time (s)');
        end
        if mod(j - 1, nCols) == 0
            ylabel('ΔR/R (%)');
        end
        box off;
    end

    sgtitle(sprintf('Sessions by Date for %s', mouse_id));
end
