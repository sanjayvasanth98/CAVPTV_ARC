function trackingResults = particle_tracking_CAVPTV(bufferedCentroids, blobMasks, parameters, folderSaveDir, fileName, inputType, videoFile, prompt, bgImage, maskArray)
    numFrames = length(bufferedCentroids);
    tracks = [];
    trackID = 1;

    % Kalman Filter parameters
    dt = 1;
    A = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];
    H = [1 0 0 0; 0 1 0 0];
    Q_default = parameters.processNoiseCovariance * eye(4);
    R = parameters.measurementNoiseCovariance * eye(2);

    % Precompute wall Y positions for each X
    [height, width] = size(maskArray);
    wallYPositions = NaN(1, width);
    for x = 1:width
        columnMask = maskArray(:, x);
        wallY = find(columnMask, 1, 'first'); % first location on x where y==1 is the wall layer
        if ~isempty(wallY)
            wallYPositions(x) = wallY;
        end
    end

    for frameIdx = 1:numFrames
        centroids = bufferedCentroids{frameIdx};
        % blobMasks(:,:,frameIdx)= blobMasks{frameIdx};
        if isempty(centroids), continue; end

        % Predict existing tracks
        for t = 1:length(tracks)
            if isPointNearWall(tracks(t).state(1:2), wallYPositions, parameters.layersNearWall)
                Q = parameters.processNoiseCovarianceNearWall * eye(4);
            else
                Q = Q_default;
            end
            [tracks(t).state, tracks(t).covariance] = kalmanPredict(tracks(t).state, tracks(t).covariance, A, Q);
        end

        % Data association
        measurements = centroids;
        costMatrix = createCostMatrix(tracks, measurements, H);

        [assignments, unassignedTracks, unassignedDetections] = assignDetectionsToTracks(costMatrix, parameters.assignmentThreshold);

        % Update tracks with assignments
        tracks = updateAssignedTracks(tracks, assignments, measurements, frameIdx, H, R, blobMasks(:,:,frameIdx), parameters);

        % Update unassigned tracks
        tracks = updateUnassignedTracks(tracks, unassignedTracks, frameIdx, parameters);

        % Create new tracks
        [tracks, trackID] = createNewTracks(tracks, unassignedDetections, measurements, frameIdx, blobMasks(:,:,frameIdx), trackID, parameters);

        % Remove terminated tracks
        tracks = tracks(~[tracks.terminated]);
    end

    % Apply filters and interpolate missing points
    tracks = applyFilters(tracks, parameters);

    % Prepare tracking results
    trackingResults = compileTrackingResults(tracks, blobMasks, maskArray);
end


%% Helper Functions

function [predictedState, predictedCovariance] = kalmanPredict(state, covariance, A, Q)
    predictedState = A * state;
    predictedCovariance = A * covariance * A' + Q;
end

function [updatedState, updatedCovariance] = kalmanUpdate(predictedState, predictedCovariance, measurement, H, R)
    K = predictedCovariance * H' / (H * predictedCovariance * H' + R);
    updatedState = predictedState + K * (measurement - H * predictedState);
    updatedCovariance = (eye(size(K,1)) - K * H) * predictedCovariance;
end

function costMatrix = createCostMatrix(tracks, measurements, H)
    numTracks = length(tracks);
    numDetections = size(measurements, 1);
    costMatrix = zeros(numTracks, numDetections);

    for i = 1:numTracks
        predictedMeasurement = H * tracks(i).state;
        for j = 1:numDetections
            costMatrix(i,j) = norm(predictedMeasurement(1:2) - measurements(j,1:2)');
        end
    end
end

function [assignments, unassignedTracks, unassignedDetections] = assignDetectionsToTracks(costMatrix, costThreshold)
    numTracks = size(costMatrix, 1);
    numDetections = size(costMatrix, 2);

    assignments = [];
    unassignedTracks = 1:numTracks;
    unassignedDetections = 1:numDetections;

    for i = 1:numTracks
        [minCost, minIdx] = min(costMatrix(i, :));
        if minCost < costThreshold
            assignments = [assignments; i, minIdx];
            unassignedDetections(unassignedDetections == minIdx) = [];
            unassignedTracks(unassignedTracks == i) = [];
        end
    end
end

function tracks = updateAssignedTracks(tracks, assignments, measurements, frameIdx, H, R, blobMask, parameters)
    for i = 1:size(assignments, 1)
        trackIdx = assignments(i,1);
        detectionIdx = assignments(i,2);
        measurement = measurements(detectionIdx, :)';

        % Check if the measurement is inside a blob mask
        if parameters.cutOffAtBlobMasks
            if isPointInBlobMask(measurement(1:2), blobMask)
                tracks(trackIdx).terminated = true;
                continue;
            end
        end

        % Kalman Update
        [updatedState, updatedCovariance] = kalmanUpdate(tracks(trackIdx).state, tracks(trackIdx).covariance, measurement, H, R);

        % Update track
        tracks(trackIdx).state = updatedState;
        tracks(trackIdx).covariance = updatedCovariance;
        tracks(trackIdx).age = tracks(trackIdx).age + 1;
        tracks(trackIdx).totalVisibleCount = tracks(trackIdx).totalVisibleCount + 1;
        tracks(trackIdx).consecutiveInvisibleCount = 0;
        tracks(trackIdx).history = [tracks(trackIdx).history; frameIdx, measurement(1)', measurement(2)'];
    end
end

function tracks = updateUnassignedTracks(tracks, unassignedTracks, frameIdx, parameters)
    for i = 1:length(unassignedTracks)
        ind = unassignedTracks(i);
        tracks(ind).age = tracks(ind).age + 1;
        tracks(ind).consecutiveInvisibleCount = tracks(ind).consecutiveInvisibleCount + 1;

        if tracks(ind).consecutiveInvisibleCount >= parameters.maxTrackAge
            tracks(ind).terminated = true;
        end
    end
end

function [tracks, newTrackID] = createNewTracks(tracks, unassignedDetections, measurements, frameIdx, blobMask, trackID, parameters)
    newTrackID = trackID;
    for i = 1:length(unassignedDetections)
        detectionIdx = unassignedDetections(i);
        measurement = measurements(detectionIdx, :)';

        % Check if the detection is inside a blob mask
        if parameters.cutOffAtBlobMasks
            if isPointInBlobMask(measurement(1:2), blobMask)
                continue;
            end
        end

        % Initialize state
        initialState = [measurement(1); measurement(2); 0; 0];
        initialCovariance = 100 * eye(4);

        % Create new track
        newTrack = struct('id', newTrackID, 'state', initialState, 'covariance', initialCovariance, ...
                          'age', 1, 'totalVisibleCount', 1, 'consecutiveInvisibleCount', 0, ...
                          'terminated', false, 'history', [frameIdx, measurement(1)', measurement(2)']);
        tracks = [tracks; newTrack];
        newTrackID = newTrackID + 1;
    end
    %disp(['Created new track with ID: ', num2str(newTrackID - 1)]);
end

function inBlob = isPointInBlobMask(point, blobMask)
    x = round(point(1));
    y = round(point(2));
    [height, width] = size(blobMask);
    if x >=1 && x <= width && y >=1 && y <= height
        inBlob = blobMask(y,x);
    else
        inBlob = false;
    end
end

function nearWall = isPointNearWall(point, wallYPositions, layersNearWall)
    x = round(point(1));
    y = round(point(2));
    [height, width] = size(wallYPositions);
    if x >=1 && x <= width && ~isnan(wallYPositions(x))
        wallY = wallYPositions(x);
        if y >= wallY - layersNearWall && y <= wallY
            nearWall = true;
        else
            nearWall = false;
        end
    else
        nearWall = false;
    end
end


function tracks = applyFilters(tracks, parameters)
    disp(['Number of tracks before filtering: ', num2str(length(tracks))]);
    % Remove tracks near edges
    if parameters.offsetFromEdges > 0
        imageWidth = parameters.im_res(2);
        imageHeight = parameters.im_res(1);
        for i = 1:length(tracks)
            positions = tracks(i).history(:,2:3);
            if any(positions(:,1) <= parameters.offsetFromEdges) || ...
               any(positions(:,2) <= parameters.offsetFromEdges) || ...
               any(positions(:,1) >= imageWidth - parameters.offsetFromEdges) || ...
               any(positions(:,2) >= imageHeight - parameters.offsetFromEdges)
                tracks(i).terminated = true;
            end
        end
        tracks = tracks(~[tracks.terminated]);
        disp(['Number of tracks after filtering: ', num2str(length(tracks))]);
    end

    % Remove vertical or crossing tracks
    if parameters.removeVerticalOrCrossingTracks
        for i = 1:length(tracks)
            positions = tracks(i).history(:,2:3);
            dx = diff(positions(:,1));
            dy = diff(positions(:,2));
            slopes = dy ./ dx;
            if any(abs(slopes) > parameters.verticalThreshold)
                tracks(i).terminated = true;
            end
        end
        tracks = tracks(~[tracks.terminated]);
    end

    % Interpolate missing points
    if parameters.interpolateMissingPoints
        for i = 1:length(tracks)
            frames = tracks(i).history(:,1);
            positions = tracks(i).history(:,2:3);
            fullFrames = frames(1):frames(end);
            if length(fullFrames) > length(frames)
                % Interpolate missing frames
                interpolatedX = interp1(frames, positions(:,1), fullFrames, 'linear');
                interpolatedY = interp1(frames, positions(:,2), fullFrames, 'linear');
                tracks(i).history = [fullFrames', interpolatedX', interpolatedY'];
            end
        end
    end
end

function trackingResults = compileTrackingResults(tracks, blobMasks, maskArray)
    trackingResults = struct('x', [], 'y', [], 'trackID', [], 'frame', [], 'spotID', [], 'uniqueTrackID', [], 'blobMasks', blobMasks, 'maskArray', maskArray);
    spotID = 1;
    for i = 1:length(tracks)
        history = tracks(i).history;
        numPoints = size(history, 1);
        trackingResults.x = [trackingResults.x; history(:,2)];
        trackingResults.y = [trackingResults.y; history(:,3)];
        trackingResults.trackID = [trackingResults.trackID; tracks(i).id * ones(numPoints,1)];
        trackingResults.frame = [trackingResults.frame; history(:,1)];
        trackingResults.spotID = [trackingResults.spotID; (spotID:spotID+numPoints-1)'];
        trackingResults.uniqueTrackID = [trackingResults.uniqueTrackID; i * ones(numPoints,1)];
        spotID = spotID + numPoints;
    end
end
