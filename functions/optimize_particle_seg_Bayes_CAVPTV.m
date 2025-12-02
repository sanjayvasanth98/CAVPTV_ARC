function fval = optimize_particle_seg_Bayes_CAVPTV(params, image, groundTruthCentroids, maxDiameter, matchingThreshold)
    % Extract parameters
    cannyLowerThreshold = params.cannyLowerThreshold;
    cannyUpperThreshold = params.cannyUpperThreshold;
    morphologicalSeFactor = params.morphologicalSeFactor;
    eccentricityThreshold = params.eccentricityThreshold;
    
    % Combine parameters into a vector
    paramVector = [cannyLowerThreshold, cannyUpperThreshold, morphologicalSeFactor, eccentricityThreshold];
    
    % Perform particle segmentation
    [~, ~, segmentedCentroids] = detect_particles(paramVector, image, maxDiameter);
    
    % Ensure centroids are 2D matrices
    if isempty(segmentedCentroids)
        segmentedCentroids = [NaN, NaN];
    elseif size(segmentedCentroids, 2) ~= 2
        segmentedCentroids = [NaN, NaN];
    end
    
    % Calculate F1 score
    [~, ~, f1] = calculate_f1_score(groundTruthCentroids, segmentedCentroids, matchingThreshold);
    f1 = min(f1, 1);
    
    % Return negative F1 score as bayesopt minimizes
    fval = -f1;
end