function particles = particle_analysis_CAVPTV(parameters, folderSaveDir, fileName, inputType, videoFile, preprocessedFrames)
       %---------------------------------------------------------------------------%            
                        % Particle analysis function for CAVPTV %
        % Input Arguments:
        %   parameters: struct containing necessary parameters
        %   folderSaveDir: directory to save results
        %   fileName: name of the video or image file (without extension)
        %   inputType: 'Images' or 'Video'
        %   videoFile: path to video file or folder containing images
        %   preprocessedFrames: if not empty, contains all frames as a 3D array

        % Outputs:
        %   particles: struct containing detected centroids, number of particles, optimized parameters
       %---------------------------------------------------------------------------%

% Extract parameters
maxDiameter = parameters.maxDiameter;
useExistingGT = logical(parameters.useExistingGT);  % Boolean input, Is ground truth data available?
groundTruthData = load(parameters.groundTruthData);  % If useExistingGT == true, file location
optimized = logical(parameters.optimized);  % Boolean input
optimizationMethod = parameters.optimizationMethod;  % 1 for GA, 2 for Bayesian
optimizedParams = parameters.optimizedParams; % If optimized == true
bufferSize = parameters.bufferSize;  % For createBufferedCentroids
matchingThreshold = parameters.matchingThreshold;  % For calculate_f1_score
displayFigures = logical(parameters.displayFigures);  % Boolean to control displaying figures
Bayesiterations = parameters.Bayesiterations; % max iterations for Bayesian optimization
GA_maxgenerations = parameters.GA_maxgenerations; % max generations for GA optimization
initialParams = parameters.initialParams;  % [cannyLowerThreshold, cannyUpperThreshold, morphologicalSeFactor, eccentricityThreshold]
fileExtension = parameters.fileExtension;  % E.g., '*.tif'
GTframenumber = parameters.groundTruthframenumber; % Ground truth frame number
binaryimgsave_old = logical(parameters.binaryimgsave_old); %Boolean input, if you want to save binary images

% Initialize outputs
particles = struct();

% Set up figure visibility
if displayFigures
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');
end

% Initialize variables to store results
detectedCentroids = {};
numParticles = [];

% Handle preprocessedFrames
if ~isempty(preprocessedFrames)
    % preprocessedFrames is not empty, use frames from it
    if ndims(preprocessedFrames) == 3
        numFrames = size(preprocessedFrames, 3); % for 2d images
    elseif ndims(preprocessedFrames) == 4
        numFrames = size(preprocessedFrames, 4); % for 3d images
    else
        error('preprocessedFrames has unexpected dimensions.');
    end
else
    % Handle inputType
    if strcmpi(inputType, 'Images')
        % Assuming videoFile is a folder containing images
        folderPath = videoFile;
        fileList = dir(fullfile(folderPath, fileExtension));
        numFrames = length(fileList);
    elseif strcmpi(inputType, 'Video')
        % videoFile is the path to the video file
        videoObj = VideoReader(videoFile);
        numFrames = floor(videoObj.Duration * videoObj.FrameRate);
        % numFrames = parameters.numof_frames;
    else
        error('Invalid inputType. Must be ''Images'' or ''Video''.');
    end
end

% Read the first frame for optimization and ground truth
if ~isempty(preprocessedFrames)
    % Use the first frame from preprocessedFrames
    if ndims(preprocessedFrames) == 3
        A = preprocessedFrames(:,:,GTframenumber);
    elseif ndims(preprocessedFrames) == 4
        A = preprocessedFrames(:,:,:,GTframenumber);
    else
        error('preprocessedFrames has unexpected dimensions.');
    end
elseif strcmpi(inputType, 'Images')
    % Read the first image
    A = imread(fullfile(folderPath, fileList(1).name));
elseif strcmpi(inputType, 'Video')
    % Read the first frame from the video
    A = read(videoObj, GTframenumber);
end

% Optimization
if optimized
    % Use the provided optimized parameters
    bestParams = optimizedParams;
    disp('Using provided optimized parameters.');
else
    % Proceed with ground truth data
    if useExistingGT
        % Use provided ground truth centroids
        centroids = groundTruthData.centroids;
        disp('Using provided ground truth centroids.');
    else
        error('No ground truth data provided. Please provide ground truth centroids in parameters.groundTruthData.');
    end

    % Display the number of ground truth particles
    numGroundTruthParticles = size(centroids, 1);
    disp(['Number of ground truth particles: ', num2str(numGroundTruthParticles)]);

    % Add a buffer around each centroid
    bufferedGroundTruth = createBufferedCentroids(centroids, bufferSize);

    % Proceed with optimization
    disp('Proceeding with parameter optimization...');

    % Run particle detection using initial parameters
    [~, ~, initialSegmentedCentroids, numInitialParticles] = detect_particles(initialParams, A, maxDiameter); % maxdiameter-specified in parameters, A- frame that you would consider for optimization

    % Ensure `initialSegmentedCentroids` has two columns (even if empty)
    if isempty(initialSegmentedCentroids)
        initialSegmentedCentroids = [NaN, NaN];
    elseif size(initialSegmentedCentroids, 2) ~= 2
        initialSegmentedCentroids = [NaN, NaN];
    end

    % Display number of detected particles before optimization
    disp(['Number of detected particles before optimization: ', num2str(numInitialParticles)]);

    % Calculate precision, recall, and F1 score before optimization
    [precisionBefore, recallBefore, f1Before] = calculate_f1_score(bufferedGroundTruth, initialSegmentedCentroids, matchingThreshold);

    % Display the precision, recall, and F1 score before optimization
    disp('Initial Detection Results (Before Optimization):');
    disp(['Precision: ', num2str(precisionBefore)]);
    disp(['Recall: ', num2str(recallBefore)]);
    disp(['F1 Score: ', num2str(f1Before)]);

    % Choose Optimization Method
    if optimizationMethod == 1
        
        % Run Genetic Algorithm (GA) optimization using F1 score as the objective function
        lb = parameters.lb;   % Lower bounds for parameters
        ub = parameters.ub;   % Upper bounds for parameters

        % Run GA
        options = optimoptions('ga', 'MaxGenerations', GA_maxgenerations, 'Display', 'iter');
        bestParams = ga(@(params) optimize_particle_seg_GA_CAVPTV(params, A, bufferedGroundTruth, maxDiameter, matchingThreshold), 4, [], [], [], [], lb, ub, [], options);
    
    elseif optimizationMethod == 2
        
        % Run Bayesian Optimization
        optVars = [
            optimizableVariable('cannyLowerThreshold', parameters.cannyLowerThresholdRange)
            optimizableVariable('cannyUpperThreshold', parameters.cannyUpperThresholdRange)
            optimizableVariable('morphologicalSeFactor', parameters.morphologicalSeFactorRange)
            optimizableVariable('eccentricityThreshold', parameters.eccentricityThresholdRange)
            ];

        results = bayesopt(@(params) optimize_particle_seg_Bayes_CAVPTV(params, A, bufferedGroundTruth, maxDiameter, matchingThreshold), ...
            optVars, 'MaxObjectiveEvaluations', Bayesiterations, 'IsObjectiveDeterministic', true, ...
            'AcquisitionFunctionName', 'expected-improvement-plus', 'Verbose', 1);

        % Get best parameters from Bayesian Optimization
        bestParams = [results.XAtMinObjective.cannyLowerThreshold, ...
            results.XAtMinObjective.cannyUpperThreshold, ...
            results.XAtMinObjective.morphologicalSeFactor, ...
            results.XAtMinObjective.eccentricityThreshold];
    else
        error('Invalid optimization method. Must be 1 (GA) or 2 (Bayesian).');
    end

    % Display final parameter values after optimization
    disp('Final Parameter Values After Optimization:');
    fprintf('Canny Lower Threshold: %.4f\n', bestParams(1));
    fprintf('Canny Upper Threshold: %.4f\n', bestParams(2));
    fprintf('Morphological Structuring Element Factor: %.2f\n', bestParams(3));
    fprintf('Eccentricity Threshold: %.2f\n', bestParams(4));
end

% Plotting results (optional)
if ~optimized
    if displayFigures
        % Plot Ground Truth Centroids
        figure;
        imshow(A, []);
        hold on;
        plot(centroids(:,1), centroids(:,2), 'go', 'MarkerSize', 5, 'LineWidth', 2);
        title('Ground Truth Centroids');

        % Plot Particles Detected Before Optimization
        figure;
        imshow(A, []);
        hold on;
        plot(initialSegmentedCentroids(:,1), initialSegmentedCentroids(:,2), 'ro', 'MarkerSize', 5, 'LineWidth', 2);
        title('Particles Detected Before Optimization');

        % Plot Particles Detected After Optimization
        [~, finalDetection, finalSegmentedCentroids, numFinalParticles] = detect_particles(bestParams, A, maxDiameter);
        figure;
        imshow(A, []);
        hold on;
        plot(finalSegmentedCentroids(:,1), finalSegmentedCentroids(:,2), 'bo', 'MarkerSize', 5, 'LineWidth', 2);
        title('Particles Detected After Optimization');
    end
end
% Process all frames
for i = 1:numFrames
    % Read each frame
    if ~isempty(preprocessedFrames)
        % Use frame from preprocessedFrames
        if ndims(preprocessedFrames) == 3
            A = preprocessedFrames(:,:,i);
        elseif ndims(preprocessedFrames) == 4
            A = preprocessedFrames(:,:,:,i);
        else
            error('preprocessedFrames has unexpected dimensions.');
        end
    elseif strcmpi(inputType, 'images')
        A = imread(fullfile(folderPath, fileList(i).name));
    elseif strcmpi(inputType, 'video')
        A = read(videoObj, i);
    end

    % Detect particles using optimized parameters
    [~, finalDetection, segmentedCentroids, numDetectedParticles] = detect_particles(bestParams, A, maxDiameter);

    % Store results
    detectedCentroids{i} = segmentedCentroids;
    numParticles(i) = numDetectedParticles;

    % Optionally, save the binary image
    if binaryimgsave_old
        % Convert detected particles into binary image (white particles on black background)
        binaryImage = imbinarize(rgb2gray(finalDetection));  % Adjust as necessary
              
        % Save the binary image
        binaryFileName = fullfile(folderSaveDir, sprintf('%s_frame_%04d_binary.tif', fileName, i));
        imwrite(binaryImage, binaryFileName);
    end    

    % Optionally, display detected particles
    if displayFigures
        figure;
        imshow(finalDetection, []);
        title(['Detected Particles in Frame ', num2str(i)]);
    end
    
    if (i<10) || (i>numFrames-10) || (mod(i, 20) == 0)        
         % Display detection summary
        disp(['Frame ', num2str(i), ': Detected ', num2str(numDetectedParticles), ' particles.']);
    end

end


disp('All frames processed and binary images saved.');

% Outputs
particles.detectedCentroids = detectedCentroids;
particles.numParticles = numParticles;
particles.optimizedParams = bestParams;

end
