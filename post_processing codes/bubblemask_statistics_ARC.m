%% analyze_bubble_masks_cluster.m
% Cluster-ready analysis of bubbleMask_all from multiple cases
% - Reads .mat files containing:
%     bubbleMask_all (Ny x Nx x Nt, logical)
%     imgHeight, imgWidth, numFrames, numFramesToProcess
%     pixels_per_m, time_between_frames
% - Computes per-case stats and saves:
%     * <caseLabel>_bubbleStats.mat
%     * Figures (PNG) into saveDir
%
% Run on cluster with something like:
%   matlab -nodisplay -nosplash -nodesktop -batch "analyze_bubble_masks_cluster"

%% ==================== USER SETTINGS ====================

saveDir = '/home/USERNAME/Inception_analysis/plots';  % <-- EDIT

cases = {
    struct('label','Smooth',  'file','/home/USERNAME/Inception_analysis/Videos/smooth/smooth_bubbleMasks.mat')
    struct('label','ks5um',   'file','/home/USERNAME/Inception_analysis/Videos/ks5/somefile_bubbleMasks.mat')
    struct('label','ks20um',  'file','/home/USERNAME/Inception_analysis/Videos/ks20/somefile_bubbleMasks.mat')
    % Add more as needed...
};

% Maximum number of frames to process per case (to limit cost on cluster)
maxFramesPerCase = Inf;   % set e.g. 5000 if needed

%% ==================== GENERAL SETUP ====================

if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% Cluster environment: no on-screen figures
set(0,'DefaultFigureVisible','off');
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesLineWidth',1.2);

numCases  = numel(cases);
caseColors = lines(numCases);   % good standard distinct colors

% For overlay plots across cases
allCase_voidFrac    = cell(numCases,1);
allCase_bubbleCount = cell(numCases,1);
allCase_diameters   = cell(numCases,1);   % equivalent diameters [m]
allCase_aspect      = cell(numCases,1);
allCase_labels      = cell(numCases,1);
allCase_dt          = zeros(numCases,1);

%% ==================== MAIN CASE LOOP ====================

for c = 1:numCases
    label   = cases{c}.label;
    matFile = cases{c}.file;
    allCase_labels{c} = label;

    fprintf('=== Case %d/%d: %s ===\n', c, numCases, label);
    fprintf('  File: %s\n', matFile);

    if ~isfile(matFile)
        warning('File not found, skipping: %s', matFile);
        continue;
    end

    % Load meta data ONLY
    meta = load(matFile, 'imgHeight','imgWidth', ...
                        'numFrames','numFramesToProcess', ...
                        'pixels_per_m','time_between_frames');
    imgHeight          = meta.imgHeight;
    imgWidth           = meta.imgWidth;
    numFrames          = meta.numFrames;
    numFramesToProcess = meta.numFramesToProcess;
    pixels_per_m       = meta.pixels_per_m;
    dt                 = meta.time_between_frames;

    allCase_dt(c) = dt;

    % Access bubbleMask_all via matfile to avoid loading entire 3D array at once
    m = matfile(matFile);
    info = whos(m, 'bubbleMask_all');
    if isempty(info)
        warning('No variable bubbleMask_all in %s. Skipping.', matFile);
        continue;
    end

    nFramesInMask = info.size(3);
    nFrames = min([nFramesInMask, numFramesToProcess, maxFramesPerCase]);
    fprintf('  Using %d frames (of %d available)\n', nFrames, nFramesInMask);

    % Preallocate frame-wise stats
    voidFrac    = zeros(nFrames,1);
    bubbleCount = zeros(nFrames,1);

    % For distributions
    allDiameters   = [];   % equivalent diameter [m]
    allAspectRatio = [];   % MajorAxisLength / MinorAxisLength

    % For intermittency map (time-averaged bubble presence)
    sumMask = zeros(imgHeight, imgWidth, 'double');

    %% -------- Frame loop for this case --------
    for k = 1:nFrames
        % Slice bubbleMask for frame k
        mask = m.bubbleMask_all(:,:,k);
        if ~islogical(mask)
            mask = logical(mask);
        end

        % Void fraction for this frame
        voidFrac(k) = nnz(mask) / numel(mask);

        % Regionprops for bubble-level geometry
        stats = regionprops(mask, 'Area','MajorAxisLength','MinorAxisLength');

        bubbleCount(k) = numel(stats);

        if ~isempty(stats)
            areasPx = [stats.Area]';   % pixels^2
            % Equivalent diameter in pixels (circle with same area)
            equivDiam_px = 2*sqrt(areasPx/pi);   % px

            % Convert to meters
            equivDiam_m = equivDiam_px / pixels_per_m;

            % Aspect ratio: Major / Minor
            maj = [stats.MajorAxisLength]';
            minax = [stats.MinorAxisLength]';
            % Avoid divide by zero
            minax(minax == 0) = eps;
            aspectRatio = maj ./ minax;

            allDiameters   = [allDiameters;   equivDiam_m];
            allAspectRatio = [allAspectRatio; aspectRatio];
        end

        % Accumulate for intermittency
        sumMask = sumMask + double(mask);
    end

    %% -------- Derive intermittency / profiles --------
    intermittency = sumMask / nFrames;  % probability of bubble presence at each pixel

    % Wall-normal void-fraction profile (average over x and time)
    alpha_y = sum(sumMask,2) / (nFrames * imgWidth);   % size: [imgHeight x 1]

    % Save stats for this case
    statsFile = fullfile(saveDir, sprintf('%s_bubbleStats.mat', label));
    save(statsFile, ...
        'voidFrac','bubbleCount', ...
        'allDiameters','allAspectRatio', ...
        'intermittency','alpha_y', ...
        'imgHeight','imgWidth','nFrames', ...
        'pixels_per_m','dt');

    fprintf('  Saved stats: %s\n', statsFile);

    % Store for cross-case plots
    allCase_voidFrac{c}    = voidFrac;
    allCase_bubbleCount{c} = bubbleCount;
    allCase_diameters{c}   = allDiameters;
    allCase_aspect{c}      = allAspectRatio;

    %% -------- Plots specific to this case --------

    % 1) Intermittency 2D map
    fig1 = figure;
    imagesc(intermittency);
    axis image;
    colormap(fig1, 'parula');
    colorbar;
    xlabel('x (pixels)');
    ylabel('y (pixels)');
    title(sprintf('Intermittency map: %s', label), 'Interpreter','none');
    set(gca,'YDir','normal','Box','on');

    print(fig1, fullfile(saveDir, sprintf('%s_intermittency.png', label)), '-dpng','-r300');
    close(fig1);

    % 2) Wall-normal void-fraction profile
    % Convert y index to meters if you want; here just pixels
    y_pix = (0:imgHeight-1)';  % pixel index from top; adjust if needed

    fig2 = figure;
    plot(alpha_y, y_pix, 'LineWidth',1.8);
    set(gca,'YDir','reverse'); % assuming y=0 at top, wall near top; adjust as needed
    xlabel('Time-averaged void fraction');
    ylabel('y (pixels)');
    title(sprintf('Void-fraction profile: %s', label), 'Interpreter','none');
    box on;

    print(fig2, fullfile(saveDir, sprintf('%s_voidFractionProfile.png', label)), '-dpng','-r300');
    close(fig2);

end % case loop

%% ==================== CROSS-CASE PLOTS (PDFs & TIME SERIES) ====================

% Helper for legend labels only for cases that actually ran
validCases = find(~cellfun(@isempty, allCase_voidFrac));
nValid = numel(validCases);
if nValid == 0
    warning('No valid cases processed. Exiting.');
    return;
end

%% 1) Void fraction time series (overlay)
figV = figure;
hold on;
for idx = 1:nValid
    cID = validCases(idx);
    t = (0:numel(allCase_voidFrac{cID})-1)' * allCase_dt(cID);
    plot(t, allCase_voidFrac{cID}, 'Color', caseColors(cID,:), 'DisplayName', allCase_labels{cID});
end
hold off;
xlabel('Time (s)');
ylabel('Void fraction (-)');
title('Global void fraction vs time');
box on;
legend('Location','best','Interpreter','none');

print(figV, fullfile(saveDir, 'allCases_voidFraction_timeSeries.png'), '-dpng','-r300');
close(figV);

%% 2) PDF of bubble count per frame
figN = figure;
hold on;
for idx = 1:nValid
    cID = validCases(idx);
    N = allCase_bubbleCount{cID};
    if isempty(N), continue; end
    maxN = max(N);
    edges = 0:1:maxN;
    if numel(edges) < 2, continue; end
    counts = histcounts(N, edges, 'Normalization','pdf');
    centers = 0.5*(edges(1:end-1) + edges(2:end));
    plot(centers, counts, '-', 'Color', caseColors(cID,:), 'DisplayName', allCase_labels{cID});
end
hold off;
xlabel('Bubble count per frame');
ylabel('PDF');
title('PDF of bubble count per frame');
box on;
legend('Location','best','Interpreter','none');

print(figN, fullfile(saveDir, 'allCases_pdf_bubbleCount.png'), '-dpng','-r300');
close(figN);

%% 3) PDF of bubble size (equivalent diameter in meters)
figD = figure;
hold on;
for idx = 1:nValid
    cID = validCases(idx);
    D = allCase_diameters{cID};
    if isempty(D), continue; end
    maxD = prctile(D, 99);   % avoid long tails from a few large bubbles
    edges = linspace(0, maxD, 40);
    counts = histcounts(D, edges, 'Normalization','pdf');
    centers = 0.5*(edges(1:end-1) + edges(2:end));
    plot(centers, counts, '-', 'Color', caseColors(cID,:), 'DisplayName', allCase_labels{cID});
end
hold off;
xlabel('Equivalent diameter (m)');
ylabel('PDF');
title('PDF of bubble size (equivalent diameter)');
box on;
legend('Location','best','Interpreter','none');

print(figD, fullfile(saveDir, 'allCases_pdf_bubbleSize.png'), '-dpng','-r300');
close(figD);

%% 4) PDF of aspect ratio
figA = figure;
hold on;
for idx = 1:nValid
    cID = validCases(idx);
    A = allCase_aspect{cID};
    if isempty(A), continue; end
    % Limit to avoid ridiculous outliers
    maxA = prctile(A, 99);
    minA = 1; % aspect ratio >= 1
    edges = linspace(minA, maxA, 40);
    counts = histcounts(A, edges, 'Normalization','pdf');
    centers = 0.5*(edges(1:end-1) + edges(2:end));
    plot(centers, counts, '-', 'Color', caseColors(cID,:), 'DisplayName', allCase_labels{cID});
end
hold off;
xlabel('Aspect ratio (Major/Minor)');
ylabel('PDF');
title('PDF of bubble aspect ratio');
box on;
legend('Location','best','Interpreter','none');

print(figA, fullfile(saveDir, 'allCases_pdf_aspectRatio.png'), '-dpng','-r300');
close(figA);

fprintf('All analysis complete. Figures and stats saved in:\n  %s\n', saveDir);
