function bufferedCentroids = createBufferedCentroids(centroids, bufferSize)

    % Ensure centroids is not empty and has two columns
    if isempty(centroids)
        disp('No centroids provided.');
        bufferedCentroids = [];
        return;
    elseif size(centroids, 2) < 2
        error('Centroids cell array must have at least two columns (X and Y coordinates).');
    end

    % Create a buffer around each centroid
    bufferedCentroids = [];
    for i = 1:size(centroids, 1)
        cx = centroids(i, 1);  % Extract X coordinate from the cell
        cy = centroids(i, 2);  % Extract Y coordinate from the cell
        [X, Y] = meshgrid(cx - bufferSize:cx + bufferSize, cy - bufferSize:cy + bufferSize);
        bufferedCentroids = [bufferedCentroids; X(:), Y(:)];
    end

end
