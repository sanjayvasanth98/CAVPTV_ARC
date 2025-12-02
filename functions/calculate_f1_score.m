function [precision, recall, f1] = calculate_f1_score(groundTruthCentroids, detectedCentroids, matchingThreshold)
    
    % Calculate distances between ground truth and detected centroids
    distances = pdist2(groundTruthCentroids, detectedCentroids);
    
    % Determine matches based on the matching threshold
    [minDistances, ~] = min(distances, [], 2);
    truePositives = sum(minDistances <= matchingThreshold);
    
    % Precision and Recall
    precision = truePositives / size(detectedCentroids, 1);
    recall = truePositives / size(groundTruthCentroids, 1);
    
    % F1 Score
    if precision + recall == 0
        f1 = 0;
    else
        f1 = 2 * (precision * recall) / (precision + recall);
    end
    
end