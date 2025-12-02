function Buffer_Xray = Buffer_XRay_CAVPTV(parameters, folderSaveDir, fileName, inputType, videoFile, preprocessedFrames, particles, Single_point_maxima)
%---------------------------------------------------------------------------%
%          Buffer X-ray Detection Function for CAVPTV
%
% Input Arguments:
%   parameters: struct containing necessary parameters, including:
%       - bufferSize: buffer level in pixels
%       - displayFigures: boolean to control figure display
%       - saveVideo: boolean to control saving the output video
%       - fileExtension: file extension for images (e.g., '*.tif')
%   folderSaveDir: directory to save results
%   fileName: name of the video or image file (without extension)
%   inputType: 'Images' or 'Video'
%   videoFile: path to video file or folder containing images
%   preprocessedFrames: if not empty, contains all frames as a 3D or 4D array
%   particles: struct containing detected centroids from particle analysis
%   Single_point_maxima: struct containing detected centroids from single point maxima
%
% Outputs:
%   Buffer_Xray: struct containing detected centroids and number of particles
%       detectedCentroids_BX: cell array of detected centroids per frame
%       numParticles_BX: array of number of particles per frame
%---------------------------------------------------------------------------%

% Extract parameters
bufferSize = parameters.bufferSize; % Buffer level in pixels
displayFigures = logical(parameters.displayFigures); % Control displaying figures
saveVideo = logical(parameters.saveVideo); % Control saving the output video
fileExtension = parameters.fileExtension; % E.g., '*.tif'

% Initialize outputs
Buffer_Xray = struct();

% Set up figure visibility
if displayFigures
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');
end

% Initialize variables to store results
detectedCentroids_BX = {};
numParticles_BX = [];

% Handle preprocessedFrames
if ~isempty(preprocessedFrames)
    if ndims(preprocessedFrames) == 3
        numFrames = size(preprocessedFrames, 3);
        [imgHeight, imgWidth, ~] = size(preprocessedFrames);
    elseif ndims(preprocessedFrames) == 4
        numFrames = size(preprocessedFrames, 4);
        [imgHeight, imgWidth, ~, ~] = size(preprocessedFrames);
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
        if numFrames > 0
            sampleFrame = imread(fullfile(folderPath, fileList(1).name));
            [imgHeight, imgWidth, ~] = size(sampleFrame);
        else
            error('No images found in the specified folder.');
        end
    elseif strcmpi(inputType, 'Video')
        % videoFile is the path to the video file
        videoObj = VideoReader(videoFile);
        numFrames = floor(videoObj.Duration * videoObj.FrameRate);
        sampleFrame = read(videoObj, 1);
        [imgHeight, imgWidth, ~] = size(sampleFrame);
    else
        error('Invalid inputType. Must be ''Images'' or ''Video''.');
    end
end

% Preallocate the output video if saving
if saveVideo
    outputVideoFile = fullfile(folderSaveDir, [fileName, '_BufferXray.avi']);
    outputVideo = VideoWriter(outputVideoFile, 'Grayscale AVI');
    % Set the frame rate (optional, can set it to match the input video)
    if exist('videoObj', 'var')
        outputVideo.FrameRate = videoObj.FrameRate;
    else
        outputVideo.FrameRate = 30; % Default frame rate
    end
    open(outputVideo);
end

% Process all frames
for i = 1:numFrames
    % Read each frame
    if ~isempty(preprocessedFrames)
        % Use frame from preprocessedFrames
        if ndims(preprocessedFrames) == 3
            frame = preprocessedFrames(:,:,i);
        elseif ndims(preprocessedFrames) == 4
            frame = preprocessedFrames(:,:,:,i);
        end
    elseif strcmpi(inputType, 'Images')
        frame = imread(fullfile(folderPath, fileList(i).name));
    elseif strcmpi(inputType, 'Video')
        frame = read(videoObj, i);
    end
    
    % Initialize black background image
    blackFrame = zeros(imgHeight, imgWidth, 'uint8');
    
    % Get centroids from particles and Single_point_maxima for this frame
    particlesCentroids = particles.detectedCentroids{i}; % Nx2 array
    spmCentroids = Single_point_maxima.detectedCentroids_SP{i}; % Mx2 array
    
    % Create buffered areas around particles centroids
    if ~isempty(particlesCentroids)
        bufferedAreas = createBufferedCentroids(particlesCentroids, bufferSize, imgHeight, imgWidth);
    else
        bufferedAreas = false(imgHeight, imgWidth);
    end
    
    % Initialize list to hold true particles
    trueParticles = [];
    
    % For each centroid in Single_point_maxima, check if it lies within any buffered area
    if ~isempty(spmCentroids)
        for j = 1:size(spmCentroids, 1)
            x = round(spmCentroids(j,1));
            y = round(spmCentroids(j,2));
            % Ensure indices are within image bounds
            x = max(1, min(x, imgWidth));
            y = max(1, min(y, imgHeight));
            if bufferedAreas(y, x)
                trueParticles = [trueParticles; x, y];
                % Mark the point on blackFrame
                blackFrame(y, x) = 255;
            end
        end
    end
    
    % Store results
    detectedCentroids_BX{i} = trueParticles;
    numParticles_BX(i) = size(trueParticles, 1);
    
    % Optionally display figures
    if displayFigures
        figure;
        imshow(blackFrame);
        title(['Buffer X-ray Points - Frame ', num2str(i)]);
    end
    
    % Write frame to video if saving
    if saveVideo
        writeVideo(outputVideo, blackFrame);
    end
end

% Close video writer if needed
if saveVideo
    close(outputVideo);
    disp(['Buffer X-ray video saved to ', outputVideoFile]);
end

% Outputs
Buffer_Xray.bufferedCentroids  = detectedCentroids_BX;
Buffer_Xray.numParticles_BX = numParticles_BX;

end

% --- Helper Function ---
function bufferedAreas = createBufferedCentroids(centroids, bufferSize, imgHeight, imgWidth)
% Create a binary mask of buffered areas around centroids
% centroids: Nx2 array of centroids [x, y]
% bufferSize: buffer radius in pixels
% imgHeight, imgWidth: dimensions of the image
% bufferedAreas: binary mask of size [imgHeight, imgWidth]

% Initialize the binary mask
bufferedAreas = false(imgHeight, imgWidth);

% For each centroid, create a circular buffer area
for i = 1:size(centroids, 1)
    xCenter = centroids(i,1);
    yCenter = centroids(i,2);
    % Create a grid of coordinates
    [X, Y] = meshgrid(1:imgWidth, 1:imgHeight);
    % Calculate distance from the centroid
    distance = sqrt((X - xCenter).^2 + (Y - yCenter).^2);
    % Mark pixels within the bufferSize
    bufferedAreas(distance <= bufferSize) = true;
end
end
