function score = optimize_particle_seg_GA_CAVPTV(params, image, groundTruthCentroids, maxDiameter, matchingThreshold)
    % Perform particle segmentation
    [~, ~, segmentedCentroids] = detect_particles(params, image, maxDiameter);
    
    % Ensure centroids are 2D matrices
    if isempty(segmentedCentroids)
        segmentedCentroids = [NaN, NaN];
    elseif size(segmentedCentroids, 2) ~= 2
        segmentedCentroids = [NaN, NaN];
    end
    
    % Calculate F1 score
    [~, ~, f1] = calculate_f1_score(groundTruthCentroids, segmentedCentroids, matchingThreshold);
    f1 = min(f1, 1);
    
    % Return negative F1 score as GA minimizes
    score = -f1;
end