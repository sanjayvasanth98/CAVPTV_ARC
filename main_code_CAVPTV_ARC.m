% ---------------------------------------------------------------------
% Main script for cavitation inception PTV analysis
% Cluster-ready (no GUI, no interactive prompts)
% Author: Sanjay Vasanth
% ---------------------------------------------------------------------

clearvars; clc; close all;

%% --------------------------------------------------------------------
% 1. USER-SET PATHS (EDIT THESE ONLY ONCE)
% --------------------------------------------------------------------

%uncomment bottom lines for adding path to toolbox on ARC

% % --- Load your personal toolbox path safely ---
% userPathFile = fullfile(getenv('HOME'), 'matlab', 'pathdef.m');
% if isfile(userPathFile)
%     addpath(genpath(fileparts(userPathFile)));  % add folder containing pathdef.m
%     run(userPathFile);                          % execute your saved pathdef.m
%     fprintf('Loaded custom MATLAB path: %s\n', userPathFile);
% else
%     warning('Custom pathdef.m not found. Using default MATLAB path.');
% end


% Folder containing all your CAVPTV functions
functionFolder = 'G:\My Drive\Research\cavitation\Final codes for analysis\CAVPTV ARC\functions_CAVPTV';    % <<< EDIT THIS >>>

% Parameters file absolute path (a .m file)
paramsFile    = 'G:\My Drive\Research\cavitation\Final codes for analysis\CAVPTV ARC\functions_CAVPTV\parameters_CAVPTV.m';  % <<< EDIT THIS >>>

% Directory where your input videos are stored
videoFolder   = 'G:\My Drive\Research\cavitation\Final codes for analysis\CAVPTV ARC\Examples\No particle inception data\New folder';          % <<< EDIT THIS >>>

% Output folder to save all results
baseSaveDir   = 'G:\My Drive\Research\cavitation\Final codes for analysis\CAVPTV ARC\Examples\No particle inception data\New folder\test';              % <<< EDIT THIS >>>

%BACKGROUND OPTION (NO PROMPTS)
% 1 = Provide a background stack and compute average
% 2 = Load single averaged background image
% 3 = No background (recommended default)
inputType = 'Video';
backgroundMode = 3;


%% --------------------------------------------------------------------
% 2. ADD FUNCTION FOLDER TO MATLAB PATH
% --------------------------------------------------------------------
addpath(functionFolder);
fprintf('Added function folder: %s\n', functionFolder);

%% --------------------------------------------------------------------
% 3. START PARALLEL POOL (cluster-safe)
% --------------------------------------------------------------------
disp('Initializing parallel pool...');
if isempty(gcp('nocreate'))
    parpool('local');   % works on SLURM/OpenMP/cluster nodes
    disp('Parallel pool initialized.');
else
    disp('Parallel pool already exists.');
end


%% --------------------------------------------------------------------
% 4. LOAD PARAMETERS (as function)
% --------------------------------------------------------------------
[~, paramFunctionName, ~] = fileparts(paramsFile);
fprintf('Using parameters function: %s\n', paramFunctionName);


%% --------------------------------------------------------------------
% 5. GET VIDEO LIST FROM FOLDER (cluster-compatible)
% --------------------------------------------------------------------
videoList = dir(fullfile(videoFolder, '*.mp4'));
if isempty(videoList)
    videoList = dir(fullfile(videoFolder, '*.avi'));
end
if isempty(videoList)
    error('No video files found in: %s', videoFolder);
end

videoFiles = fullfile({videoList.folder}, {videoList.name})';
fprintf('Found %d videos.\n', numel(videoFiles));


%% --------------------------------------------------------------------
% 6. LOAD BACKGROUND IMAGES (no GUI, auto-selection)
% --------------------------------------------------------------------
background_image_extension = '*.tif';
bgImages = selectBackgroundImage_CAVPTV(inputType, videoFiles, background_image_extension, backgroundMode);



%% --------------------------------------------------------------------
% 7. CREATE MASK ARRAY
% --------------------------------------------------------------------
maskArray = maskarray_CAVPTV(videoFiles, baseSaveDir, bgImages);


%% --------------------------------------------------------------------
% 8. READ ANALYSIS TYPE AND SAMPLE TYPE FROM PARAMETERS
% --------------------------------------------------------------------
% Load default parameters from first video
parameters = feval(paramFunctionName, videoFiles{1});

analysisType = parameters.analysisType;   % must be defined in param file
sampleType   = parameters.sampleType;     % must be defined in param file


%% --------------------------------------------------------------------
% 9. BATCH PROCESSING LOOP FOR ALL VIDEOS
% --------------------------------------------------------------------
for video_num = 1:length(videoFiles)

    videoFile = videoFiles{video_num};

    % Reload parameters for each video
    parameters = feval(paramFunctionName, videoFile);

    % Make folder for this video
    [~, fileName] = fileparts(videoFile);
    folderSaveDir = fullfile(baseSaveDir, fileName);

    if ~exist(folderSaveDir, 'dir')
        mkdir(folderSaveDir);
    end

    fprintf('\n-------------------------------------------------------\n');
    fprintf('Processing video %d/%d: %s\n', video_num, length(videoFiles), videoFile);
    fprintf('Saving results to: %s\n', folderSaveDir);
    fprintf('-------------------------------------------------------\n');

    try
        tic
        if sampleType == 1 || sampleType == 2
            [detectedCentroids_SP, Single_point_maxima, bubbleMask_all, ...
             blobMasks, bufferedcentroids] = ...
                execute_analysis_CAVPTV( ...
                    sampleType, analysisType, parameters, ...
                    folderSaveDir, fileName, 'Video', videoFile, ...
                    0, bgImages{video_num}, maskArray(:,:,video_num));
        else
            execute_analysis_CAVPTV( ...
                    sampleType, analysisType, parameters, ...
                    folderSaveDir, fileName, 'Video', videoFile, ...
                    0, bgImages{video_num}, maskArray(:,:,video_num), ...
                    video_num);
        end
        toc

    catch ME
        fprintf('Error while processing %s:\n  %s\n', videoFile, ME.message);
        rethrow(ME);
    end

end

fprintf('\nBatch processing completed for all videos.\n');
