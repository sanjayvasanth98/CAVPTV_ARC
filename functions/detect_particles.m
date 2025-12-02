function [originalImage, finalDetection, segmentedCentroids, numDetectedParticles] = detect_particles(params, image, maxDiameter)
    % Extract parameters
    cannyLowerThreshold = params(1);
    cannyUpperThreshold = params(2);
    morphologicalSeFactor = params(3);
    eccentricityThreshold = params(4);
    
    % Convert the input image to grayscale if it's RGB
    if size(image, 3) == 3
        X = rgb2gray(image);
    else
        X = image;
    end
    X = im2double(X);
    
    % Edge Detection and Morphological Filtering
    I2 = edge(X, 'Canny', [cannyLowerThreshold cannyUpperThreshold]);
    I3 = imfill(I2, 'holes');
    se = strel('disk', round(maxDiameter / morphologicalSeFactor));
    I4 = imopen(I3, se);
    I51 = bwpropfilt(I4, 'MajorAxisLength', [0 maxDiameter]);
    I61 = bwpropfilt(I51, 'Eccentricity', [0 eccentricityThreshold]);
    
    % Watershed Segmentation for Overlapping Particles
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
        segmentedCentroids = [NaN, NaN];
    elseif size(segmentedCentroids, 2) ~= 2
        segmentedCentroids = [NaN, NaN];
    end
    
    % Get the number of detected particles
    numDetectedParticles = size(segmentedCentroids, 1);
    
    % Overlay results for final detection display
    originalImage = image;
    finalDetection = imoverlay(image, bwperim(finalParticles), [0 1 0]);
end