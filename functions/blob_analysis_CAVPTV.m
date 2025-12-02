function blob_analysis_CAVPTV(parameters, folderSaveDir, fileName, inputType, inputPath, preprocessedFrames, maskArray, video_num, varargin)
%------------------------------------------------------------------------%
% % Function to perform blob analysis % %
%
% This function processes a video or image sequence to compute optical flow,
% segment bubbles, calculate velocities, and save the results. In addition,
% it saves two video files:
%   (1) A binary video of the bubble masks.
%   (2) A preprocessed video (grayscale) of the input frames.
%
% Additionally, if a mask region is provided via maskArray, that mask is 
% subtracted from the bubble mask in every frame.
%
% Arguments:
%   parameters:         struct containing all the parameters (e.g., startFrame,
%                       endFrame, dt, px, verbose, etc.)
%   folderSaveDir:      Path to save processed data.
%   fileName:           Base name of the file to be saved.
%   inputType:          'Video' or 'Images'.
%   inputPath:          Path to the video file or image sequence folder.
%   preprocessedFrames: If not empty, contains all frames as a 3D (grayscale)
%                       or 4D (color) array.
%   maskArray:          A 3D logical array containing mask regions. The mask
%                       for this video is taken from maskArray(:,:,video_num).
%                       (Pass [] if not used.)
%   video_num:          The index into maskArray along the third dimension.
%   varargin:           Additional optional arguments.
%------------------------------------------------------------------------%

% ----------------- BASIC PARAMETERS ------------------------------------
startFrame          = parameters.startFrame;   % default = 1
endFrame            = parameters.endFrame;     % default = numFrames
verbose             = logical(parameters.verbose);
time_between_frames = parameters.dt;
pixels_per_m        = parameters.px;

% Flag to control optical flow / velocity computation & saving
if isfield(parameters, 'computeFlow')
    computeFlow = logical(parameters.computeFlow);
else
    computeFlow = true;  % default: ON
end

% ----------------- OPEN INPUT / DIMENSIONS -----------------------------
if ~isempty(preprocessedFrames)
    if ndims(preprocessedFrames) == 3
        [imgHeight, imgWidth, numFrames] = size(preprocessedFrames);
    elseif ndims(preprocessedFrames) == 4
        [imgHeight, imgWidth, ~, numFrames] = size(preprocessedFrames);
    else
        error('preprocessedFrames has unexpected dimensions.');
    end
else
    if strcmpi(inputType, 'Video')
        videoFile = inputPath;
        vidReader = VideoReader(videoFile);
        numFrames = vidReader.NumFrames;
        frame = readFrame(vidReader);
        frameGray = im2double(im2gray(frame));
        [imgHeight, imgWidth] = size(frameGray);
        vidReader.CurrentTime = 0;
    elseif strcmpi(inputType, 'Images')
        folder = inputPath;
        imageFiles = dir(fullfile(folder, '*.tif'));
        numFrames = length(imageFiles);
        frame = imread(fullfile(folder, imageFiles(1).name));
        frameGray = im2double(im2gray(frame));
        [imgHeight, imgWidth] = size(frameGray);
    else
        error('Invalid input type. Must be ''Video'' or ''Images''.');
    end
end

% Adjust frame indices if necessary
if isempty(endFrame) || endFrame > numFrames
    endFrame = numFrames;
end
numFramesToProcess = endFrame - startFrame + 1;

% ----------------- PREALLOCATIONS --------------------------------------

% Always store bubble masks and preprocessed frames
bubbleMask_all        = false(imgHeight, imgWidth, numFramesToProcess);
preprocessedVideo_all = zeros(imgHeight, imgWidth, numFramesToProcess, 'uint8');

% Only allocate flow-related arrays if computeFlow is true
if computeFlow
    u_real_all             = zeros(imgHeight, imgWidth, numFramesToProcess, 'single');
    v_real_all             = zeros(imgHeight, imgWidth, numFramesToProcess, 'single');
    velocity_magnitude_all = zeros(imgHeight, imgWidth, numFramesToProcess, 'single');
    globalMinVelocity      = inf;
    globalMaxVelocity      = -inf;
else
    u_real_all             = [];
    v_real_all             = [];
    velocity_magnitude_all = [];
    globalMinVelocity      = [];
    globalMaxVelocity      = [];
end

% Initialize optical flow object only if needed
if computeFlow
    opticFlow = opticalFlowRAFT;
else
    opticFlow = [];
end

if verbose
    disp('Running Blob Analysis...');
    if computeFlow
        disp('  Optical flow: ON');
    else
        disp('  Optical flow: OFF');
    end
end

% ----------------- MAIN LOOP -------------------------------------------
frameCounter = 1;
for frameIndex = startFrame:endFrame

    % --- Read / pick frame ---
    if ~isempty(preprocessedFrames)
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

    % store preprocessed grayscale
    preprocessedVideo_all(:,:,frameCounter) = uint8(255 * frameGray);

    % --- Segment bubbles ---
    bubbleMask = segmentBubbles_CAVPTV(frameGray, parameters);

    % --- Apply maskArray region (if any) ---
    if ~isempty(maskArray)
        if iscell(maskArray)
            maskRegion = maskArray{video_num};
        else
            maskRegion = maskArray(:,:,video_num);
        end

        if ~isempty(maskRegion)
            if ~islogical(maskRegion)
                maskRegion = logical(maskRegion);
            end
            bubbleMask = bubbleMask & ~maskRegion;
        end
    end

    % --- Optical flow & velocities (optional) ---
    if computeFlow
        flow = estimateFlow(opticFlow, frameGray);

        uMasked = flow.Vx .* bubbleMask;
        vMasked = flow.Vy .* bubbleMask;

        u_real = (uMasked / pixels_per_m) / time_between_frames;
        v_real = (vMasked / pixels_per_m) / time_between_frames;
        velocity_magnitude = sqrt(u_real.^2 + v_real.^2);

        velocity_magnitude_bubbles = velocity_magnitude(bubbleMask);
        if ~isempty(velocity_magnitude_bubbles)
            globalMinVelocity = min(globalMinVelocity, min(velocity_magnitude_bubbles(:)));
            globalMaxVelocity = max(globalMaxVelocity, max(velocity_magnitude_bubbles(:)));
        end

        % store flow results
        u_real_all(:,:,frameCounter)             = single(u_real);
        v_real_all(:,:,frameCounter)             = single(v_real);
        velocity_magnitude_all(:,:,frameCounter) = single(velocity_magnitude);
    end

    % store mask
    bubbleMask_all(:,:,frameCounter) = bubbleMask;

    frameCounter = frameCounter + 1;
end

% ----------------- SAVE RESULTS ----------------------------------------

% 1) Save flow results (only if computeFlow == true)
if computeFlow
    saveFileName = fullfile(folderSaveDir, strcat(fileName, '_flow.mat'));
    save(saveFileName, 'u_real_all', 'v_real_all', 'velocity_magnitude_all', ...
        'bubbleMask_all', 'globalMinVelocity', 'globalMaxVelocity', ...
        'imgHeight', 'imgWidth', 'numFrames', 'numFramesToProcess', ...
        'pixels_per_m', 'time_between_frames', '-v7.3');
end

% 2) Save bubble masks + calibration always (no flow required)
bubbleMaskMatName = fullfile(folderSaveDir, strcat(fileName, '_bubbleMasks.mat'));
save(bubbleMaskMatName, 'bubbleMask_all', ...
    'imgHeight', 'imgWidth', 'numFrames', 'numFramesToProcess', ...
    'pixels_per_m', 'time_between_frames', '-v7.3');


% 3) Save Bubble Mask Video (Grayscale)
aviFileName = fullfile(folderSaveDir, strcat(fileName, '_bubbleMasks.avi'));
writerObj = VideoWriter(aviFileName);
writerObj.FrameRate = 1;  % adjust if needed
open(writerObj);
for i = 1:numFramesToProcess
    maskFrame = uint8(bubbleMask_all(:,:,i)) * 255;
    writeVideo(writerObj, maskFrame);
end
close(writerObj);

% 4) Save Preprocessed Video (Grayscale)
preprocAviFileName = fullfile(folderSaveDir, strcat(fileName, '_preprocessed.avi'));
writerObjPreproc = VideoWriter(preprocAviFileName);
writerObjPreproc.FrameRate = 1;
open(writerObjPreproc);
for i = 1:numFramesToProcess
    preprocFrame = preprocessedVideo_all(:,:,i);
    writeVideo(writerObjPreproc, preprocFrame);
end
close(writerObjPreproc);

if verbose
    fprintf('Blob analysis ended. Results saved in: %s\n', folderSaveDir);
    if computeFlow
        fprintf('  Flow results: %s\n', saveFileName);
    else
        fprintf('  Flow results: SKIPPED (computeFlow = false)\n');
    end
    fprintf('  Bubble masks MAT: %s\n', bubbleMaskMatName);
    fprintf('  Bubble masks AVI: %s\n', aviFileName);
    fprintf('  Preprocessed AVI: %s\n', preprocAviFileName);
end

end
