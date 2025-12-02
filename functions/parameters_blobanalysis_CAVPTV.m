function parameters = parameters_CAVPTV(videoFile,optimizedParams)
    % Define the parameters structure with all relevant settings

    %% Video Parameters
    parameters.acq_mode = 'time_resolved'; % image acquisition mode % not used

    % Dynamically determine camera resolution from the video file
    videoObj = VideoReader(videoFile);
    parameters.im_res = [videoObj.Height, videoObj.Width]; % Camera resolution from the video
    parameters.numof_frames = videoObj.NumFrames;

    parameters.n_frames = 9; % number of double-frames 
    % parameters.dt = 3.68e-6; % time separation of image pairs [s]
    % parameters.m = 0.0015625; % mapping scale [mm/px]
    % parameters.px = 640000; % mapping scale [px/m]
    % parameters.im_pre = 'I'; % image file prefix
    % parameters.im_file = 'tif'; % image file format
    %load('im_roi.mat'); 
    %parameters.im_roi = im_roi; % import image mask, if defined
    % parameters.im_roi = ones(parameters.im_res(1), parameters.im_res(2)); % define mask else

    %% File Selection Parameters
    parameters.analysis_mode = 2;  % 1 for single file, 2 for batch processing
    parameters.analysis_type = 1;  % 1 for Particle Analysis, 2 for Blob Analysis
    parameters.single_file = 'path/to/single/file.avi';  % Single file path
    parameters.batch_files = {'path/to/file1.avi', 'path/to/file2.avi'};  % Batch files
    parameters.save_dir = 'path/to/save/directory';  % Save directory

    %% Processing/Tracking Parameters
    % Particle detection
    % parameters.p_size = 3; % Particle size in px
    % parameters.p_int  = 2; % Particle intensity

   %% Defining analysis and sampletype

   % Define Analysis Type
    parameters.analysisType = 1; % '1' for Particle Analysis, '2' for Blob Analysis 
   % Define sampleType
    parameters.sampleType = 1;  % '1' for X-ray Samples, '2' for Backlight Samples with particles, '3' for Backlight Samples Without Particles

    %% Background image
    % parameters.image_extension = '*.tif'; % Specify the desired image extension

    %% Pre-processing 

    parameters.median_filter_size = [10 10]; % Size of the median filter
    parameters.image_extension = '*.tif';    % File extension for images when inputType is 'Images'
    parameters.lowerbound = 110; % Lower percentile for contrast stretching
    parameters.higherbound = 145;% Upper percentile for contrast stretching
    parameters.savevideopreprocess = false;
    parameters.xrayillumination = false; % logical input for indicating that the type of illumination for process variations
    
    % Preprocessing Control
    parameters.preprocessingEnabled = true; % true to enable preprocessing, false to disable
    
    %% particle analysis
    parameters.maxDiameter = 10;  % Set the maximum equivalent diameter
    parameters.useExistingGT = true;  % Set to true if ground truth data is available
    parameters.groundTruthData = "C:\Users\kbsanjayvasanth\OneDrive - Virginia Tech\Research\cavitation\PTV\Particle detection and tracking Sanjay\Final\ground_truth_centroids.mat";  % Provide the ground truth centroids
    parameters.optimized = true;  % Set to true if you have optimized parameters
    parameters.optimizationMethod = 2;  % 1 for GA, 2 for Bayesian
    parameters.bufferSize = 3;  % Adjust as necessary
    parameters.matchingThreshold = 5;  % Adjust as necessary
    parameters.displayFigures = true;  % Set to true to display figures
    parameters.fileExtension = '*.tif';  % File extension for images
    parameters.initialParams = [0.01, 0.05, 12, 0.85];  % Initial parameters for optimization
    parameters.lb = [0.001, 0.05, 3, 0.5];  % Lower bounds for GA -obtained from unit testing
    parameters.ub = [0.01, 0.15, 15, 0.9];  % Upper bounds for GA -obtained from unit testing
    parameters.cannyLowerThresholdRange = [0.001, 0.01];  % For Bayesian optimization -obtained from unit testing
    parameters.cannyUpperThresholdRange = [0.05, 0.15]; % For Bayesian optimization -obtained from unit testing
    parameters.morphologicalSeFactorRange = [3, 15]; % For Bayesian optimization -obtained from unit testing
    parameters.eccentricityThresholdRange = [0.5, 0.9]; % For Bayesian optimization -obtained from unit testing
    parameters.Bayesiterations = 50;   % max iterations for Bayesian optimization
    parameters.GA_maxgenerations = 50; % max generations for GA optimization
    parameters.groundTruthframenumber = 1; %ground truth frame number in the whole video eg: frame 1- 1 ,frame 2- 2
    parameters.binaryimgsave_old = false; %vBoolean input, if you want to save binary images

    % Check if optimized parameters are provided
    if nargin == 2  % If optimizedParams is passed
        parameters.optimizedParams = optimizedParams;
    else
        parameters.optimizedParams = [0.009135913764441,0.080614869392425,8.760271147312201,0.832773865256616];  % Parameters after optimization
    end

    %% Single point maxima detection - X-ray Illumination (comment out this section if testing Backlight illumination samples)

    parameters.tolerance = 98; %tolerance for finding the particle
    parameters.strict = false; % Use strict mode
    parameters.lightbackground = false; % Background type, Light - true, Dark - false
    parameters.threshold = -Inf; % Minimum value for maxima to be considered (-Inf for no threshold)
    parameters.outputType = 'SINGLE_POINTS'; % Output format
    parameters.excludeOnEdges = true; % Exclude detection at the edges
    parameters.fileExtension = '*.tif'; % File extension for images
    parameters.saveMaximaVideo = true; % Save the video? true/false

    %% Single point maxima detection - Backlight Illumination (comment out this section if testing X-ray illumination samples)

    % parameters.tolerance = 5; %tolerance for finding the particle
    % parameters.strict = false; % Use strict mode
    % parameters.lightbackground = true; % Background type, Light - true, Dark - false
    % parameters.threshold = Inf; % Minimum value for maxima to be considered (-Inf for no threshold)
    % parameters.outputType = 'SINGLE_POINTS'; % Output format
    % parameters.excludeOnEdges = true; % Exclude detection at the edges
    % parameters.fileExtension = '*.tif'; % File extension for images
    % parameters.saveMaximaVideo = true; % Save the video? true/false

    % buffer
    parameters.combinemasks = true; % do you want to combine masks of Venturi and blobs?

    %% X-Ray buffer

    parameters.bufferSize = 3; % Buffer size in pixels
    parameters.saveVideo = true;
    parameters.fileExtension = '*.tif';

    %% Tracking Parameters
    parameters.processNoiseCovariance = 50;              % Scalar value for process noise covariance Q
    parameters.measurementNoiseCovariance = 50;          % Scalar value for measurement noise covariance R
    parameters.assignmentThreshold = 50;                % Maximum acceptable cost for assignment
    parameters.maxTrackAge = 2;                         % Maximum frames a track can be inactive
    parameters.cutOffAtBlobMasks = true;                % Cut off tracks at blob masks
    parameters.offsetFromEdges = 5;                    % Offset distance from edges in pixels
    parameters.removeVerticalOrCrossingTracks = true;   % Remove vertical or crossing tracks
    parameters.verticalThreshold = 2;                   % Threshold for vertical movement
    parameters.interpolateMissingPoints = true;         % Interpolate missing points in tracks
    parameters.layersNearWall = 65; % Number of layers near the wall
    parameters.processNoiseCovarianceNearWall = 1; %eg: 1e-2 % Adjusted process noise near the wall

    %% blob analysis
    parameters.startFrame =1; %modify this to start processing frames from this number
    parameters.endFrame = [] ; %modify this to end processing frames at this number
    parameters.verbose = true;
    % parameters.dt = 3.68e-6; % Time between frames in seconds (X-ray)
    % parameters.px = 640000;  % Pixels per meter (X-ray)
    parameters.dt = 3.71e-6; % Pixels per meter (Backlight)
    parameters.px = 256000; % Pixels per meter (Backlight)
    parameters.xrayill = 3; % 1=Xray images 2=backlight images with particles 3=backlight images without particles 4= segementing bubbles from backlight videos, 5 = cloud cavitation
    %% Double-Frame Processing
    parameters.track_method = 'hist_match'; % not used
    parameters.f_o_s = 35; % not used
    parameters.n_neighbours = 5; % not used
    parameters.gauss_interp = 1; % not used

    %% Time-resolved Processing
    parameters.min_dist = 5; % not used
    parameters.n_mp_ti = 1; % not used

    %% Additional Parameters
    parameters.sigma_threshold = 10; % not used

    % Save the parameters structure dynamically after setting im_res
    save('parameters_CAVPTV.mat', 'parameters');
end