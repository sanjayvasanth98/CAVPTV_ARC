function [bubbleMask_all,blobMasks]=blob_analysismaskreturn_CAVPTV(parameters, folderSaveDir, fileName, inputType, inputPath, preprocessedFrames, varargin)
%------------------------------------------------------------------------%
                % % Function to perform blob analysis % %
    % This function processes a video or image sequence to compute 
    % optical flow, segment bubbles, calculate velocities, and 
    % save the results.
    
    % Arguments:
    %   parameters: struct containing all the parameters
    %   folderSaveDir: Path to save processed data
    %   fileName: name of the file
    %   inputType: 'video' or 'images'
    %   inputPath: Path to the video file or image sequence folder
    %   preprocessedFrames: if not empty, contains all frames as a 3D or 4D array

%------------------------------------------------------------------------%

startFrame          = parameters.startFrame; % 'StartFrame' (default = 1): The first frame to process.
endFrame            = parameters.endFrame;   % 'EndFrame' (default = numFrames): The last frame to process.
verbose             = logical(parameters.verbose);     % 'Verbose' (default = true): Display processing messages.
time_between_frames = parameters.dt;
pixels_per_m        = parameters.px;

% Initialize variables to store global min and max velocity magnitudes
globalMinVelocity = inf;
globalMaxVelocity = -inf;

% Handle preprocessedFrames
if ~isempty(preprocessedFrames)
    % preprocessedFrames is not empty, use frames from it
    if ndims(preprocessedFrames) == 3
        % Grayscale frames
        [imgHeight, imgWidth, numFrames] = size(preprocessedFrames);
    elseif ndims(preprocessedFrames) == 4
        % Color frames
        [imgHeight, imgWidth, ~, numFrames] = size(preprocessedFrames);
    else
        error('preprocessedFrames has unexpected dimensions.');
    end
else
    % Determine input type and set up accordingly
    if strcmpi(inputType, 'Video')
        % Video processing
        videoFile = inputPath;
        vidReader = VideoReader(videoFile);
        numFrames = vidReader.NumFrames;

        % Read the first frame to get image size
        frame = readFrame(vidReader);
        frameGray = im2double(im2gray(frame));
        [imgHeight, imgWidth] = size(frameGray);

        % Reset the reader to the first frame
        vidReader.CurrentTime = 0;
    elseif strcmpi(inputType, 'Images')
        % Image sequence processing
        folder = inputPath;
        imageFiles = dir(fullfile(folder, '*.tif'));
        numFrames = length(imageFiles);

        % Read the first image to get image size
        frame = imread(fullfile(folder, imageFiles(1).name));
        frameGray = im2double(im2gray(frame));
        [imgHeight, imgWidth] = size(frameGray);
    else
        error('Invalid input type. Must be ''Video'' or ''Images''.');
    end
end

% Adjust the number of frames and frame indices
if isempty(endFrame) || endFrame > numFrames
    endFrame = numFrames;
end
numFramesToProcess = endFrame - startFrame + 1;

% Adjust the arrays to store data
u_real_all = zeros(imgHeight, imgWidth, numFramesToProcess, 'single');
v_real_all = zeros(imgHeight, imgWidth, numFramesToProcess, 'single');
velocity_magnitude_all = zeros(imgHeight, imgWidth, numFramesToProcess, 'single');
bubbleMask_all = false(imgHeight, imgWidth, numFramesToProcess);

% Initialize the opticFlow object
opticFlow = opticalFlowRAFT;
blobMasks = cell(numFramesToProcess, 1); % Updated to Cell Array

if verbose
    disp('Running Blob Analysis...');
end

% Process each frame
frameCounter = 1;
for frameIndex = startFrame:endFrame
    % Read frame
    if ~isempty(preprocessedFrames)
        % Use frame from preprocessedFrames
        if ndims(preprocessedFrames) == 3
            frame = preprocessedFrames(:,:,frameIndex);
        elseif ndims(preprocessedFrames) == 4
            frame = preprocessedFrames(:,:,:,frameIndex);
        else
            error('preprocessedFrames has unexpected dimensions.');
        end
        frameGray = im2double(im2gray(frame));
    elseif strcmpi(inputType, 'Video')
        vidReader.CurrentTime = (frameIndex - 1) / vidReader.FrameRate;
        frame = readFrame(vidReader);
        frameGray = im2double(im2gray(frame));
    elseif strcmpi(inputType, 'Images')
        frame = imread(fullfile(folder, imageFiles(frameIndex).name));
        frameGray = im2double(im2gray(frame));
    end

    % Estimate optical flow
    flow = estimateFlow(opticFlow, frameGray);

    % Bubble segmentation
    bubbleMask = segmentBubbles_CAVPTV(frameGray,parameters);

    % Calculate velocities within bubbles
    uMasked = flow.Vx .* bubbleMask;
    vMasked = flow.Vy .* bubbleMask;

    % Convert to real-world units (m/s)
    u_real = (uMasked / pixels_per_m) / time_between_frames;
    v_real = (vMasked / pixels_per_m) / time_between_frames;

    % Calculate velocity magnitude within bubbles
    velocity_magnitude = sqrt(u_real.^2 + v_real.^2);
    velocity_magnitude_bubbles = velocity_magnitude(bubbleMask);

    % Update global min and max velocities
    if ~isempty(velocity_magnitude_bubbles)
        globalMinVelocity = min(globalMinVelocity, min(velocity_magnitude_bubbles(:)));
        globalMaxVelocity = max(globalMaxVelocity, max(velocity_magnitude_bubbles(:)));
    end

    % Store data
    u_real_all(:,:,frameCounter) = single(u_real);
    v_real_all(:,:,frameCounter) = single(v_real);
    velocity_magnitude_all(:,:,frameCounter) = single(velocity_magnitude);
    bubbleMask_all(:,:,frameCounter) = bubbleMask;
    blobMasks{frameCounter} = bubbleMask;

    frameCounter = frameCounter + 1;
end

% Save all data into a single .mat file using the passed 'fileName'
saveFileName = fullfile(folderSaveDir, strcat(fileName, '.mat'));

% Extract the folder path from the save file name
folderPath = fileparts(saveFileName);

% Check if the folder exists; if not, create it
if ~exist(folderPath, 'dir')
    mkdirStatus = mkdir(folderPath);
    if ~mkdirStatus
        error('Failed to create the directory: %s', folderPath);
    else
        fprintf('Directory created: %s\n', folderPath);
    end
end

disp(saveFileName);

% Check if the directory was successfully created
if ~exist(folderPath, 'dir')
    error('The directory still does not exist after attempting to create it: %s', folderPath);
end

if length(saveFileName) > 260
    warning('The file path exceeds 260 characters: %d characters', length(saveFileName));
end

save(saveFileName, 'u_real_all', 'v_real_all', 'velocity_magnitude_all', 'bubbleMask_all', ...
    'globalMinVelocity', 'globalMaxVelocity', 'imgHeight', 'imgWidth', 'numFrames', ...
    'pixels_per_m', 'time_between_frames','-v7.3');

if verbose
    fprintf('Blobs processing ended. Results saved in: %s\n', folderSaveDir);
end
end
