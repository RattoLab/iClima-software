%% script of the GCaMP6f spectra
%% INPUT
% ratio830nm_45=ratio;
% ratio920nm_45=ratio;
% ratio=ratio830nm_45;
% ratio=ratio920nm_45;
% Time=t;
%% %% ---- 1) Load .analyzed file-----
[fileName, pathName] = uigetfile('*analyzed*.mat', 'Select analyzed file');
if isequal(fileName,0)
    error('No file selected');
end
data = load(fullfile(pathName, fileName));% Load the .mat file
fn = fieldnames(data);
S = data.(fn{1});

% Ensure Time is column vector
Time = S.Time(:);

ratio = S.G_sub_dark;
% ratio = S.Ratio;% Assume S.Ratio is [time x trials x ROIs]

dt = mean(diff(Time));% Frame period
%% Assume ratio is time x trials x ROIs
[nTime, nTrials, nROIs] = size(ratio);

% STEP 1 — F0 from pre-stimulus time points
baselineIdx = Time < 0;

F0 = nan(1, nTrials, nROIs);
for tr = 1:nTrials
    for roi = 1:nROIs
        F0(1,tr,roi) = mean(ratio(baselineIdx, tr, roi), 'omitnan');
    end
end

% STEP 2 — ΔF/F0 (percent)
dFoverF = (ratio - F0) ./ F0;
dFoverF = dFoverF * 100;   % %

% STEP 3 — Plot each trial
mov_avg_window = [5 0];  % moving average window

if nROIs == 1
    % Special case: only one ROI
    nRows = ceil(sqrt(nTrials));
    nCols = ceil(nTrials / nRows);
    
    figure('Name', 'All Trials (Single ROI)', 'NumberTitle', 'off');
    tl = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for tr = 1:nTrials
        ax = nexttile;
        hold(ax,'on')
        
        y = movmean(dFoverF(:, tr, 1), mov_avg_window, 'omitnan');
        plot(ax, Time, y, 'k', 'LineWidth', 1.5);
        
        % Stimulus patch
        stimDur = 5;
        yl = [-10, 3];
        patch(ax, [0 stimDur stimDur 0], [yl(1) yl(1) yl(2) yl(2)], ...
            [0.45 0.55 0.45], 'EdgeColor','none', 'FaceAlpha',0.4);
        
        % Stimulus marker
        plot(ax, [0 0], yl, 'r--', 'LineWidth', 1)
        
        % Axes formatting
%         xlim(ax, [-20 25])
        ylim(ax, yl)
         xlim(ax, 'auto')
        ylim(ax, 'auto')
   %         set(ax,'ydir','reverse')
        title(ax, ['Trial ' num2str(tr)])
        box(ax,'off')
    end
    
    xlabel(tl, 'Time (s)')
    ylabel(tl, '\DeltaR/R (%)')
    
else
    % Multiple ROIs: original plotting per trial
    for tr = 1:nTrials
        nRows = ceil(sqrt(nROIs));
        nCols = ceil(nROIs / nRows);

        figure('Name', ['Trial ' num2str(tr)], 'NumberTitle', 'off');
        tl = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

        for roi = 1:nROIs
            ax = nexttile;
            hold(ax,'on')

            y = movmean(dFoverF(:,tr,roi), mov_avg_window, 'omitnan');
            plot(ax, Time, y, 'k', 'LineWidth', 1.5);

            % Stimulus patch
            stimDur = 5;
            yl = [-15, 10];
            patch(ax, [0 stimDur stimDur 0], [yl(1) yl(1) yl(2) yl(2)], ...
                [0.45 0.55 0.45], 'EdgeColor','none', 'FaceAlpha',0.4);

            % Stimulus marker
            plot(ax, [0 0], yl, 'r--', 'LineWidth', 1)

            % Axes formatting
%             xlim(ax, [-20 25])
            ylim(ax, yl)
             xlim(ax, 'auto')
             ylim(ax, 'auto')
%             set(gca,'ydir','reverse')
            title(ax, ['ROI ' num2str(roi)])
            box(ax,'off')
        end

        xlabel(tl, 'Time (s)')
        ylabel(tl, '\DeltaR/R (%)')
    end
end
%% %% Plot mean ± SEM across trials for each ROI
nROIs   = size(ratio, 3);
nTrials = size(ratio, 2);

Time = S.Time(:);           
baselineIdx = Time < 0;

nCols = ceil(sqrt(nROIs));
nRows = ceil(nROIs / nCols);
meanAllDeltaF   = nan(length(Time), nROIs);
figure('Name', 'mean ± SEM across trials for each ROI', 'NumberTitle', 'off')
tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact')
ROIs=[1,3,17,20,21,27,28]
for roi =1:nROIs
    nexttile
    hold on

    Green =squeeze(ratio(:,:,roi));   % [time × trials]
     % Baseline normalization (per trial)
    Green0 = mean(Green(baselineIdx, :), 1);
    dFoverF = 100 .* (Green - Green0) ./ Green0;

    mov_avg_window = [5 0]; 
    dFoverF = movmean(  dFoverF, mov_avg_window, 'omitnan');
    % Mean and SEM
   meanDF = mean(dFoverF, 2); 
   semDgreen  = std(dFoverF, 0, 2) ./ sqrt(nTrials);
   meanAllDeltaF(:,roi) = meanDF;
   fill([Time; flipud(Time)], ...
         [meanDF-semDgreen; flipud(meanDF+semDgreen)], [0.4 0.4 0.4], 'EdgeColor','none', 'FaceAlpha',0.4);

    % Mean trace
    plot(Time, meanDF, 'k', 'LineWidth', 1.5);
    title(sprintf('ROI %d', roi))
    xLimits = [min(Time) max(Time)];
    yLimits = [min(meanDF(:)) max(meanDF(:))];  % Y now covers all ΔRatio values
    xlim(xLimits)
    ylim( yLimits)
   xPatch = [0 5 5 0];           % fixed window
    yPatch = [yLimits(1) yLimits(1) yLimits(2) yLimits(2)];  % spans all ROIs
    patch( xPatch, yPatch, [0.9 0.6 0.9], ...
          'FaceAlpha',0.3, 'EdgeColor','none')
    if roi > (nRows-1)*nCols
        xlabel('Time')
    end
    if mod(roi-1, nCols) == 0
        ylabel('\DeltaF/F (%)')
    end
%     set(gca,'ydir','reverse')
    box off
    hold off
end
%% NEW SECTION — Filter ROIs + Global Summary Plot

stimDur = 5;                                    
stimIdx = Time >= 0 & Time <=2;

% ---- STEP 1: Remove ROIs with positive median during stimulus ----
roiKeep = false(1,nROIs);

for roi = 1:nROIs
    roiData = squeeze(dFoverF(:,:,roi));   % time x trials
    stimMedian = median(roiData(stimIdx,:), 'all','omitnan');
    if stimMedian <= 0
        roiKeep(roi) = true;
    end
end

roi_removed_idx = find(~roiKeep);
roi_kept_idx    = find(~roiKeep);

dFoverF_clean   = dFoverF(:,:,~roiKeep);
nROIs_clean     = numel(roi_kept_idx);
nROIs_removed   = numel(roi_removed_idx);

% ---- Display counts in command window ----
if isempty(roi_removed_idx)
    disp(['Kept: ' num2str(nROIs_clean) ' | Removed: 0 | Removed ROIs: none'])
else
    disp(['Kept: ' num2str(nROIs_clean) ...
          ' | Removed: ' num2str(nROIs_removed) ...
          ' | Removed ROIs: ' num2str(roi_removed_idx)])
end

% ---- STEP 2: Compute averages ----
mov_avg_window = [5 0];  % moving average window

% Average across trials (per ROI) and smooth each ROI
roi_meanTrials_raw = squeeze(mean(dFoverF_clean, 2, 'omitnan'));  % time x ROIs
roi_meanTrials = nan(size(roi_meanTrials_raw));
for roi = 1:nROIs_clean
    roi_meanTrials(:,roi) = movmean(roi_meanTrials_raw(:,roi), mov_avg_window, 'omitnan');
end

% Global mean: average of smoothed ROI traces
global_mean = mean(roi_meanTrials, 2, 'omitnan');

% Average across ROIs (per trial) for bottom subplot
trial_meanROIs = mean(dFoverF_clean, 3, 'omitnan');  % time x trials

% ---- STEP 3: Plot summary figure ----
figure('Name','Filtered ROI Summary','NumberTitle','off');
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

yl = [-6 6];
xLimits = [Time(1) Time(end)];

% -------- TOP: All ROIs (trial-averaged & smoothed) --------
ax1 = nexttile;
hold(ax1,'on')

% Plot each ROI (smoothed)
for roi = 1:nROIs_clean
    plot(ax1, Time, roi_meanTrials(:,roi), 'Color',[0.75 0.75 0.75])
end

% Plot global mean (already smoothed)
plot(ax1, Time, global_mean, 'k','LineWidth',2)

% Stimulus patch + vertical line
patch(ax1, [0 stimDur stimDur 0], [yl(1) yl(1) yl(2) yl(2)], ...
      [0.85 0.85 0.45], 'EdgeColor','none','FaceAlpha',0.3);
plot(ax1, [0 0], yl, 'r--')

% xPositions = [0, 5,10, 15,20,25,30,35];  % example x positions for vertical lines
% for x = xPositions
%     xline(ax1, x, 'r--', 'LineWidth', 1);
% end
%     xPositions = [0, 5];
%     for x = xPositions
%         if x == 25
%             xline(ax1, x, 'b--', 'LineWidth', 1.5);
%         else
%             xline(ax1, x, 'r--', 'LineWidth', 1);
%         end
%     end

% Axes formatting
xlim(ax1, xLimits)
ylim(ax1, yl)
set(ax1,'ydir','reverse')
ylabel(ax1,'\DeltaR/R (%)')

% Multi-line title: kept & removed info
% Top subplot title: first line = kept/removed info, second line = descriptive subtitle
title(ax1, {...
    ['Kept: ' num2str(nROIs_clean) ' | Removed: ' num2str(nROIs_removed) ' | Removed ROI(s): ' num2str(roi_removed_idx)], ...
    'All ROIs (Trial-averaged)'});
box(ax1,'off')

% -------- BOTTOM: All Trials (ROI-averaged) --------
ax2 = nexttile;
hold(ax2,'on')

% Optionally, smooth each trial-averaged trace
for tr = 1:nTrials
    y = movmean(trial_meanROIs(:,tr), mov_avg_window, 'omitnan');
    plot(ax2, Time, y, 'Color',[0.75 0.75 0.75])
end

% Plot global mean (from smoothed ROI traces)
plot(ax2, Time, global_mean,'k','LineWidth',3)

% Stimulus patch + vertical line
patch(ax2, [0 stimDur stimDur 0], [yl(1) yl(1) yl(2) yl(2)], ...
      [0.85 0.85 0.45], 'EdgeColor','none','FaceAlpha',0.3);
plot(ax2, [0 0], yl, 'r--')

% xPositions = [0, 5,10, 15,20,25,30,35];  % example x positions for vertical lines
% for x = xPositions
%     xline(ax2, x, 'r--', 'LineWidth', 1);
% end
%     xPositions = [0, 5,10, 15,20,25,30,35,40,45,50];
%     for x = xPositions
%         if x == 25
%             xline(ax2, x, 'b--', 'LineWidth', 1.5);
%         else
%             xline(ax2, x, 'r--', 'LineWidth', 1);
%         end
%     end

% Axes formatting
xlim(ax2,xLimits)
ylim(ax2,yl)
set(ax2,'ydir','reverse')
xlabel(ax2,'Time (s)')
ylabel(ax2,'\DeltaR/R (%)')
title(ax2,'All Trials (ROI-averaged)')
box(ax2,'off')


%% Assume ratio is time x trials x ROIs

% STEP 3 — Concatenate trials for each ROI
% Create a new "time" vector for concatenated trials
trialDur = max(Time) - min(Time);
concatTime = (0:nTime*nTrials-1)' * (trialDur/(nTime-1));  % assumes uniform spacing

% STEP 4 — Plot each ROI in separate tiles
nRows = ceil(sqrt(nROIs));
nCols = ceil(nROIs / nRows);

figure('Name','All ROIs concatenated','NumberTitle','off');
tl = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

mov_avg_window = [5 0];

for roi = 1:nROIs
    ax = nexttile;
    hold(ax,'on')
    
    % Concatenate all trials for this ROI
    yConcat = [];
    for tr = 1:nTrials
        y = movmean(dFoverF(:,tr,roi), mov_avg_window, 'omitnan');
        yConcat = [yConcat; y];  %#ok<AGROW>
    end
    
    plot(ax, concatTime, yConcat, 'k', 'LineWidth', 1.5);
     xlim(ax, [min(concatTime) max(concatTime)])
     ylim(ax, yl)
    % Optional: draw vertical lines separating trials
    for tr = 1:nTrials-1
        xSep = tr * trialDur;
        plot(ax, [xSep xSep], ylim(ax), 'r-', 'LineWidth', 1);
    end
    
    xlabel(ax,'Time (s)')
    ylabel(ax,'\DeltaR/R (%)')
    set(gca,'ydir','reverse')
    title(ax, ['ROI ' num2str(roi)])
    box(ax,'off')
end
