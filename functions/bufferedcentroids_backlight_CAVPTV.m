function bufferedcentroids_backlight = bufferedcentroids_backlight_CAVPTV(centroids, blobMasks, parameters, folderSaveDir, fileName, maskArray)

%---------------------------------------------------------------------------%
    % bufferedcentroids_backlight_CAVPTV Removes centroids lying within blob masks.
%
%   bufferedcentroids_backlight = bufferedcentroids_backlight_CAVPTV(centroids, blobmasks, parameters, folderSaveDir, fileName, maskArray)
%
%   Inputs:
%       centroids      - Cell array of detected centroids per frame (from Singlepoint_maxima_CAVPTV)
%       blobmasks      - Cell array of binary masks per frame (from blobmasks)
%       parameters     - Struct containing necessary parameters, including 'combinemasks'
%       folderSaveDir  - Directory path to save the results
%       fileName       - Base name for the saved file
%       maskArray      - Cell array of additional masks per frame (if applicable)
%
%   Output:
%       bufferedcentroids_backlight - Cell array of filtered centroids per frame
%---------------------------------------------------------------------------%
    % Validate Inputs
    if nargin < 6
        error('All six input arguments are required.');
    end

    if ~iscell(centroids)
        error('centroids must be a cell array.');
    end

    if ~iscell(blobMasks)
        error('blobMasks must be a cell array.');
    end

    %---------------------------------------------------------%
    
    combineMasks = logical(parameters.combinemasks);

    numFrames = length(centroids);
    if length(blobMasks) ~= numFrames
        error('centroids, blobmasks must have the same number of frames.');
    end

    % Initialize Output
    bufferedcentroids_backlight = cell(numFrames, 1);

    % Create Save Directory if it does not exist
    if ~exist(folderSaveDir, 'dir')
        mkdir(folderSaveDir);
    end

    % Initialize Video Writer
    videoSavePath = fullfile(folderSaveDir, [fileName, '_centroids_video.avi']);
    frameRate = 30; % Default frame rate
    if isfield(parameters, 'frameRate')
        frameRate = parameters.frameRate;
    end
    videoWriter = VideoWriter(videoSavePath);
    videoWriter.FrameRate = frameRate;
    open(videoWriter);

    % Process Each Frame
    for i = 1:numFrames
        currentCentroids = centroids{i};          % Nx2 array [x, y]
        currentBlobMask = blobMasks{i};            % Binary mask
        currentMaskArray = maskArray;           % Additional binary mask (if applicable)

        % Determine Combined Mask Based on 'combinemasks' Parameter
        if combineMasks
            if ~isempty(currentMaskArray)
                combinedMask = currentBlobMask | currentMaskArray;
            else
                combinedMask = currentBlobMask;
            end
        else
            combinedMask = currentBlobMask;
        end

        % Initialize logical array to keep track of centroids to keep
        keepCentroid = true(size(currentCentroids, 1), 1);

        % Get Image Dimensions
        [imgHeight, imgWidth] = size(combinedMask);

        % Extract all centroid coordinates and round them
        xCoords = round(currentCentroids(:,1));
        yCoords = round(currentCentroids(:,2));

        % Remove centroids that are out of bounds
        validIndices = (xCoords >=1 & xCoords <= imgWidth & yCoords >=1 & yCoords <= imgHeight);
        keepCentroid = keepCentroid & validIndices;

        % Update coordinates to only valid ones
        validX = xCoords(validIndices);
        validY = yCoords(validIndices);

        % Linear Indices of valid centroids
        linearIdx = sub2ind(size(combinedMask), validY, validX);

        % Determine which centroids lie within the combined mask
        inMask = combinedMask(linearIdx);

        % Mark centroids within the mask for removal
        keepCentroid(validIndices) = ~inMask;

        % Filter Centroids
        bufferedcentroids_backlight{i} = currentCentroids(keepCentroid, :);

        % Display Progress
        % fprintf('Processed frame %d/%d: Original Centroids = %d, Filtered Centroids = %d\n', ...
            % i, numFrames, size(currentCentroids, 1), size(bufferedcentroids_backlight{i}, 1));

        % Create a blank image (black background)
        frameImage = zeros(imgHeight, imgWidth, 3, 'uint8'); % RGB image

        % Plot centroids as white pixels
        if ~isempty(bufferedcentroids_backlight{i})
            % Ensure coordinates are within bounds
            validPlotX = bufferedcentroids_backlight{i}(:,1);
            validPlotY = bufferedcentroids_backlight{i}(:,2);

            % Round coordinates to nearest integer
            validPlotX = round(validPlotX);
            validPlotY = round(validPlotY);

            % Remove any out-of-bounds points
            validPoints = (validPlotX >=1 & validPlotX <= imgWidth & validPlotY >=1 & validPlotY <= imgHeight);
            validPlotX = validPlotX(validPoints);
            validPlotY = validPlotY(validPoints);

            % Set the pixel locations to white
            linearIndices = sub2ind([imgHeight, imgWidth], validPlotY, validPlotX);
            frameImage(linearIndices) = 255; % White color for all RGB channels
            frameImage(linearIndices + imgHeight*imgWidth) = 255;
            frameImage(linearIndices + 2*imgHeight*imgWidth) = 255;
        end

        % Write the frame to the video
        writeVideo(videoWriter, frameImage);
    end

    % Close the Video Writer
    close(videoWriter);

    % Prepare Save Filename for Centroids
    saveFileName = fullfile(folderSaveDir, [fileName, '_bufferedcentroids_backlight.mat']);

    % Save the Filtered Centroids
    save(saveFileName, 'bufferedcentroids_backlight');

    fprintf('Filtered centroids saved to %s\n', saveFileName);
    fprintf('Centroids video saved to %s\n', videoSavePath);
end
