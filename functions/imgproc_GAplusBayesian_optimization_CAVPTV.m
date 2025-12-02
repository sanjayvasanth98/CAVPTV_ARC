% Clear workspace and command window
clear all; 
clc;
close all;

% Set up some default figure and display properties
set(0, 'DefaultFigureVisible', 'on');
set(gcf, 'Color', 'none');
set(gca, 'box', 'off');

% Prompt the user for the maximum equivalent diameter of the particles
maxDiameter = input('Enter the maximum equivalent diameter of the particles to detect (in pixels): ');

% Specify the path to the folder containing the TIFF images
folderPath = 'C:\Users\kbsanjayvasanth\OneDrive - Virginia Tech\Research\cavitation\APS Data Processed\Ra5\Test\test2\testsep24';
fileExtension = '*.tif';  % File extension for TIFF images
saveDir = uigetdir([], 'Select Directory to Save Results');

% Get a list of file names matching the specified pattern
fileList = dir(fullfile(folderPath, fileExtension));

% Read the first image (used in the GA optimization and for ground truth)
A = imread(fullfile(folderPath, fileList(1).name));

% Prompt user if they have an existing ground truth file
useExistingGT = input('Do you have a ground truth file? (1 for Yes, 0 for No): ');

% Specify directory to save or upload ground truth data
if useExistingGT
    % If the user has a ground truth file, prompt for the file location
    [file, path] = uigetfile('*.mat', 'Select the Ground Truth File');
    groundTruthFile = fullfile(path, file);
    load(groundTruthFile, 'centroids');  % Load existing ground truth centroids
    disp('Using uploaded ground truth centroids.');
else
    % If no ground truth exists, prompt the user to create one
    disp('You will be prompted to select the centroids of particles.');
    
    % Show the image for manual ground truth creation
    figure;
    imshow(A, []);
    title('Click on the centroids of all particles, press Enter when done');
    hold on;
    
    % Allow the user to click the centroids and collect the points
    [centroidX, centroidY] = getpts;
    centroids = [centroidX, centroidY];
    
    % Plot the selected centroids
    plot(centroidX, centroidY, 'b+', 'MarkerSize', 10, 'LineWidth', 2);
    
    % Prompt user for the directory to save the ground truth data
    saveDir = uigetdir([], 'Select Directory to Save Ground Truth Data');
    groundTruthFile = fullfile(saveDir, 'ground_truth_centroids.mat');
    
    % Save the ground truth centroids
    save(groundTruthFile, 'centroids');
    disp('Ground truth centroids saved.');
end

% Display the number of ground truth particles
numGroundTruthParticles = size(centroids, 1);
disp(['Number of ground truth particles: ', num2str(numGroundTruthParticles)]);

% Add a buffer of a few pixels around each centroid to account for slight deviations
bufferSize = 3;  % Adjust this value based on the typical particle size
bufferedGroundTruth = createBufferedCentroids(centroids, bufferSize);

% Ask the user if parameters have already been optimized
optimized = input('Have you already optimized the parameters? (1 for Yes, 0 for No): ');

if optimized
    % If parameters are already optimized, proceed to detect particles for each frame
    disp('Proceeding with existing optimized parameters...');

    % Ask the user for optimized parameters
    cannyLowerThreshold = input('Enter the optimized Canny lower threshold: ');
    cannyUpperThreshold = input('Enter the optimized Canny upper threshold: ');
    morphologicalSeFactor = input('Enter the optimized morphological structuring element factor: ');
    eccentricityThreshold = input('Enter the optimized eccentricity threshold: ');
    bestParams = [cannyLowerThreshold, cannyUpperThreshold, morphologicalSeFactor, eccentricityThreshold];
    
else
    % If parameters have not been optimized, proceed with optimization
    disp('Proceeding with parameter optimization...');
    
    %% Step 1: Run the initial particle detection (before optimization) and compare with ground truth
    % Initial parameters
    cannyLowerThreshold = 0.01;
    cannyUpperThreshold = 0.05;
    morphologicalSeFactor = 12;
    eccentricityThreshold = 0.85;

    % Run particle detection using initial parameters
    [~, finalDetection, initialSegmentedCentroids, numInitialParticles] = detect_particles([cannyLowerThreshold, cannyUpperThreshold, morphologicalSeFactor, eccentricityThreshold], A, maxDiameter);

    % Ensure `initialSegmentedCentroids` has two columns (even if empty)
    if isempty(initialSegmentedCentroids)
        initialSegmentedCentroids = [NaN, NaN];
    elseif size(initialSegmentedCentroids, 2) ~= 2
        initialSegmentedCentroids = [NaN, NaN];
    end

    % Display number of detected particles before optimization
    disp(['Number of detected particles before optimization: ', num2str(numInitialParticles)]);

    % Calculate precision, recall, and F1 score before optimization
    [precisionBefore, recallBefore, f1Before] = calculate_f1_score(bufferedGroundTruth, initialSegmentedCentroids);

    % Display the precision, recall, and F1 score before optimization
    disp('Initial Detection Results (Before Optimization):');
    disp(['Precision: ', num2str(precisionBefore)]);
    disp(['Recall: ', num2str(recallBefore)]);
    disp(['F1 Score: ', num2str(f1Before)]);

    %% Step 2: Choose Optimization Method (GA or Bayesian)
    optimizationMethod = input('Which optimization method would you like to use? (1 for GA, 2 for Bayesian): ');

    if optimizationMethod == 1
        %% Step 2A: Run Genetic Algorithm (GA) optimization using F1 score as the objective function
        lb = [0.001, 0.05, 3, 0.5];   % Lower bounds for parameters
        ub = [0.01, 0.15, 15, 0.9];   % Upper bounds for parameters

        % Run GA
        options = optimoptions('ga', 'MaxGenerations', 50, 'Display', 'iter');
        bestParams = ga(@(params) optimize_particle_detection(params, A, bufferedGroundTruth, maxDiameter), 4, [], [], [], [], lb, ub, [], options);

    elseif optimizationMethod == 2
        %% Step 2B: Run Bayesian Optimization using the same objective function
        optVars = [
            optimizableVariable('cannyLowerThreshold', [0.001, 0.01])
            optimizableVariable('cannyUpperThreshold', [0.05, 0.15])
            optimizableVariable('morphologicalSeFactor', [3, 15])
            optimizableVariable('eccentricityThreshold', [0.5, 0.9])
        ];

        results = bayesopt(@(params) optimize_particle_detection_bayes(params, A, bufferedGroundTruth, maxDiameter), ...
            optVars, 'MaxObjectiveEvaluations', 50, 'IsObjectiveDeterministic', true, ...
            'AcquisitionFunctionName', 'expected-improvement-plus', 'Verbose', 1);

        % Get best parameters from Bayesian Optimization
        bestParams = [results.XAtMinObjective.cannyLowerThreshold, ...
                      results.XAtMinObjective.cannyUpperThreshold, ...
                      results.XAtMinObjective.morphologicalSeFactor, ...
                      results.XAtMinObjective.eccentricityThreshold];
    end

    % Display final parameter values after optimization
    disp('Final Parameter Values After Optimization:');
    fprintf('Canny Lower Threshold: %.4f\n', bestParams(1));
    fprintf('Canny Upper Threshold: %.4f\n', bestParams(2));
    fprintf('Morphological Structuring Element Factor: %.2f\n', bestParams(3));
    fprintf('Eccentricity Threshold: %.2f\n', bestParams(4));
end

%% Plot Results: Ground Truth, Before Optimization, and After Optimization

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

%% Step 3: Detect and save particles for all frames
for i = 1:length(fileList)
    % Read each image
    A = imread(fullfile(folderPath, fileList(i).name));
    
    % Detect particles using optimized parameters
    [~, finalDetection, ~, numDetectedParticles] = detect_particles(bestParams, A, maxDiameter);
    
    % Display detected particles
    figure;
    imshow(finalDetection, []);
    title(['Detected Particles in Frame ', num2str(i)]);
    
    % Convert detected particles into binary image (white particles on black background)
    binaryImage = im2bw(finalDetection);
    
    % Save the binary image
    binaryFileName = fullfile(saveDir, sprintf('frame_%04d_binary.tif', i));
    imwrite(binaryImage, binaryFileName);
    
    % Display detection summary
    disp(['Frame ', num2str(i), ': Detected ', num2str(numDetectedParticles), ' particles.']);
end

% After processing all frames
disp('All frames processed and binary images saved.');

%% Function to detect particles using the given parameters
function [originalImage, finalDetection, segmentedCentroids, numDetectedParticles] = detect_particles(params, image, maxDiameter)
    % Extract parameters
    cannyLowerThreshold = params(1);
    cannyUpperThreshold = params(2);
    morphologicalSeFactor = params(3);
    eccentricityThreshold = params(4);
    
    % Convert the input image to double format for processing
    X = im2double(image);

    %% Step 1: Edge Detection and Morphological Filtering
    I2 = edge(X, 'Canny', [cannyLowerThreshold cannyUpperThreshold]);
    I3 = imfill(I2, 'holes');
    se = strel('disk', round(maxDiameter / morphologicalSeFactor));
    I4 = imopen(I3, se);
    I51 = bwpropfilt(I4, 'MajorAxisLength', [0 maxDiameter]);
    I61 = bwpropfilt(I51, 'Eccentricity', [0 eccentricityThreshold]);

    %% Step 2: Watershed Segmentation for Overlapping Particles
    distanceTransform = bwdist(~I61);
    smoothedDistance = imgaussfilt(distanceTransform, 2);  
    L = watershed(-smoothedDistance);  
    I_watershed = I61 & (L > 0);  
    finalParticles = bwareaopen(I_watershed, 1);  

    % Get centroids of the detected particles
    stats = regionprops(finalParticles, 'Centroid');
    segmentedCentroids = cat(1, stats.Centroid);

    % Ensure segmented centroids are 2D
    if isempty(segmentedCentroids)
        segmentedCentroids = [NaN, NaN];  % Default to prevent errors if no particles detected
    elseif size(segmentedCentroids, 2) ~= 2
        segmentedCentroids = [NaN, NaN];
    end

    % Get the number of detected particles
    numDetectedParticles = size(segmentedCentroids, 1);

    % Step 4: Overlay results for final detection display
    originalImage = image;  % Original image (without overlays)
    finalDetection = imoverlay(image, bwperim(finalParticles), [0 1 0]);  % Green outline for final detection
end

%% Function to calculate precision, recall, and F1 score
function [precision, recall, f1] = calculate_f1_score(groundTruthCentroids, detectedCentroids)
    % Ensure that both groundTruthCentroids and detectedCentroids are 2D
    if size(groundTruthCentroids, 2) ~= 2
        error('Ground truth centroids must have 2 columns representing x and y.');
    end
    if size(detectedCentroids, 2) ~= 2
        error('Detected centroids must have 2 columns representing x and y.');
    end
    
    % Calculate distances between ground truth and detected centroids
    distances = pdist2(groundTruthCentroids, detectedCentroids);
    
    % Set a threshold for matching (within a few pixels)
    matchingThreshold = 5;  % Adjust as necessary
    [minDistances, ~] = min(distances, [], 2);
    
    % Precision: how many detected particles are correct
    truePositives = sum(minDistances <= matchingThreshold);
    precision = truePositives / size(detectedCentroids, 1);
    
    % Recall: how many ground truth particles were detected
    recall = truePositives / size(groundTruthCentroids, 1);
    
    % F1 Score
    if precision + recall == 0
        f1 = 0;
    else
        f1 = 2 * (precision * recall) / (precision + recall);
    end
end

%% Function to create a buffer around centroids for comparison
function bufferedCentroids = createBufferedCentroids(centroids, bufferSize)
    bufferedCentroids = [];
    for i = 1:size(centroids, 1)
        cx = centroids(i, 1);
        cy = centroids(i, 2);
        % Create a buffer around the centroid
        [X, Y] = meshgrid(cx-bufferSize:cx+bufferSize, cy-bufferSize:cy+bufferSize);
        bufferedCentroids = [bufferedCentroids; X(:), Y(:)];
    end
end

%% Function to optimize particle detection using GA
function score = optimize_particle_detection(params, image, groundTruthCentroids, maxDiameter)
    % Parameters from GA
    cannyLowerThreshold = params(1);
    cannyUpperThreshold = params(2);
    morphologicalSeFactor = params(3);
    eccentricityThreshold = params(4);
    
    % Perform particle segmentation
    [~, ~, segmentedCentroids] = detect_particles(params, image, maxDiameter);
    
    % Ensure centroids are 2D matrices for pdist2
    if isempty(segmentedCentroids)
        segmentedCentroids = [NaN, NaN];  % Ensure pdist2 doesn't fail if no particles are detected
    elseif size(segmentedCentroids, 2) ~= 2
        segmentedCentroids = [NaN, NaN];
    end
    
    % Calculate precision, recall, and F1 score
    [precision, recall, f1] = calculate_f1_score(groundTruthCentroids, segmentedCentroids);
    
    % Ensure F1 score doesn't exceed 1
    f1 = min(f1, 1);

    % Return negative F1 score as GA minimizes
    score = -f1;
end

%% Function to optimize particle detection using Bayesian optimization
function fval = optimize_particle_detection_bayes(params, image, groundTruthCentroids, maxDiameter)
    % Extract parameter values from Bayesian optimization variables
    cannyLowerThreshold = params.cannyLowerThreshold;
    cannyUpperThreshold = params.cannyUpperThreshold;
    morphologicalSeFactor = params.morphologicalSeFactor;
    eccentricityThreshold = params.eccentricityThreshold;
    
    % Combine parameters into a vector for reuse
    paramVector = [cannyLowerThreshold, cannyUpperThreshold, morphologicalSeFactor, eccentricityThreshold];
    
    % Perform particle segmentation
    [~, ~, segmentedCentroids] = detect_particles(paramVector, image, maxDiameter);
    
    % Ensure centroids are 2D matrices for pdist2
    if isempty(segmentedCentroids)
        segmentedCentroids = [NaN, NaN];  % Ensure pdist2 doesn't fail if no particles are detected
    elseif size(segmentedCentroids, 2) ~= 2
        segmentedCentroids = [NaN, NaN];
    end
    
    % Calculate precision, recall, and F1 score
    [precision, recall, f1] = calculate_f1_score(groundTruthCentroids, segmentedCentroids);
    
    % Ensure F1 score doesn't exceed 1
    f1 = min(f1, 1);

    % Return negative F1 score as bayesopt minimizes
    fval = -f1;
end
