function preprocessedFrames = preprocessing_CAVPTV(parameters, prompt, folderSaveDir, inputType, filename, videoFile, backgroundImage)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % PREPROCESSING_CAVPTV Preprocesses video frames: background
    % subtraction, median filtering, and contrast enhancement.
    %
    % Syntax:
    %   processedFrames = preprocessing_CAVPTV(parameters, prompt, folderSaveDir, inputType, filename, videoFile, backgroundImage)
    %
    % Inputs:
    %   parameters      - Structure containing various processing parameters.
    %   prompt          - Boolean flag (1 or 0) to control display of intermediate results.
    %   folderSaveDir   - Directory path where the processed outputs will be saved.
    %   inputType       - 'Images' or 'Video'.
    %   filename        - Name of the video or image file (without extension).
    %   videoFile       - Path to video file or folder containing images.
    %   backgroundImage - (Optional) Background image as a 2D matrix or path to the image file.
    %
    % Outputs:
    %   processedFrames - 3D matrix containing the processed video frames.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

displayFigures      = false; % FORCE OFF for cluster
savevideopreprocess = logical(parameters.savevideopreprocess);
xrayillumination    = logical(parameters.xrayillumination);
flatfieldcorr       = logical(parameters.flatfieldcorr);
if ~ismember(inputType, {'Images', 'Video'})
    error('Invalid inputType. Must be ''Images'' or ''Video''.');
end

if ~isfolder(folderSaveDir)
    mkdir(folderSaveDir);
    fprintf('Created save directory: %s\n', folderSaveDir);
end

%% Step 1: Read frames
if prompt
    fprintf('Reading and pre-processing frames...\n');
end

switch inputType
    case 'Images'
        folderPath   = videoFile;
        fileExtension = parameters.image_extension; % '*.png', '*.tif', etc.
        fileList = dir(fullfile(folderPath, fileExtension));
        if isempty(fileList)
            error('No images found in folder %s with extension %s.', folderPath, fileExtension);
        end
        numFrames = numel(fileList);

        firstFrame = imread(fullfile(folderPath, fileList(1).name));
        if size(firstFrame, 3) == 3
            firstFrame = im2gray(firstFrame);
        end
        [height, width] = size(firstFrame);

        videoFrames = zeros(height, width, numFrames, 'uint8');
        videoFrames(:,:,1) = firstFrame;

        for i = 2:numFrames
            frame = imread(fullfile(folderPath, fileList(i).name));
            if size(frame, 3) == 3
                frame = im2gray(frame);
            end
            videoFrames(:,:,i) = frame;
        end

    case 'Video'
        if ~isfile(videoFile)
            error('Video file does not exist: %s', videoFile);
        end
        videoObj = VideoReader(videoFile);
        numFrames = floor(videoObj.Duration * videoObj.FrameRate);
        [height, width] = deal(videoObj.Height, videoObj.Width);

        videoFrames = zeros(height, width, numFrames, 'uint8');

        frameIdx = 1;
        while hasFrame(videoObj) && frameIdx <= numFrames
            frame = readFrame(videoObj);
            grayFrame = im2gray(frame);
            videoFrames(:,:,frameIdx) = grayFrame;
            frameIdx = frameIdx + 1;
        end
end

preprocessedFrames = videoFrames;


%% Step 2: Flat-field correction (brightfield, whole stack using imflatfield)

if flatfieldcorr
    if ~isfield(parameters, 'sigma')
        error('parameters.sigma must be defined for flat-field correction.');
    end
    sigma = parameters.sigma;

    if prompt
        fprintf('Applying flat-field correction using imflatfield (sigma = %.3f)...\n', sigma);
    end

    for i = 1:numFrames
        % imflatfield works with uint8/double; preprocessedFrames is uint8 already
        preprocessedFrames(:,:,i) = imflatfield(preprocessedFrames(:,:,i), sigma);
    end
end
%% Step 3: Background image logic
% Case 3 (no BG): backgroundImage is empty → NO BG subtraction
if ~isempty(backgroundImage)
    % Background provided → use it
    if ischar(backgroundImage) || isstring(backgroundImage)
        if ~isfile(backgroundImage)
            error('Background image file does not exist: %s', backgroundImage);
        end
        avgImage = imread(backgroundImage);
        if size(avgImage, 3) == 3
            avgImage = im2gray(avgImage);
        end
    elseif isnumeric(backgroundImage)
        avgImage = backgroundImage;
        if size(avgImage, 3) == 3
            avgImage = im2gray(avgImage);
        end
    else
        error('backgroundImage must be a file path or numeric image matrix.');
    end

    if size(avgImage,1) ~= height || size(avgImage,2) ~= width
        error('Background image dimensions (%d x %d) do not match frame dimensions (%d x %d).', ...
            size(avgImage,1), size(avgImage,2), height, width);
    end

    if prompt
        fprintf('Using provided background image for subtraction.\n');
        fprintf('Performing background subtraction...\n');
    end

    % Background-subtracted stack
    k1 = 1;
    k2 = mean(avgImage(:));

    for i = 1:numFrames
        preprocessedFrames(:,:,i) = (double(preprocessedFrames(:,:,i)) - double(avgImage)) * k1 + k2;
    end

    preprocessedFrames = uint8(max(min(preprocessedFrames, 255), 0));

    % Save the background image (not a display)
    bgImagePath = fullfile(folderSaveDir, [filename '_bgimage.tif']);
    try
        imwrite(avgImage, bgImagePath);
        if prompt
            fprintf('Background image saved to %s\n', bgImagePath);
        end
    catch ME
        fprintf('Warning: Could not save background image: %s\n', ME.message);
    end
else
    % Case 3: no background image selected → NO subtraction
    if prompt
        fprintf('No background image provided → Skipping background subtraction.\n');
    end
end

%% Step 4: Median subtraction (optional, xrayillumination)
if xrayillumination
    if prompt
        fprintf('Applying median filtering and median subtraction...\n');
    end

    k1_med = 1;
    k2_med = mean(preprocessedFrames(:));

    filteredFrames = zeros(size(preprocessedFrames), 'double');
    for i = 1:numFrames
        filteredFrames(:,:,i) = medfilt2(preprocessedFrames(:,:,i), parameters.median_filter_size);
    end

    for i = 1:numFrames
        preprocessedFrames(:,:,i) = (double(preprocessedFrames(:,:,i)) - double(filteredFrames(:,:,i))) * k1_med + k2_med;
    end

    preprocessedFrames = uint8(max(min(preprocessedFrames, 255), 0));
end

%% Step 5: Contrast enhancement (optional, xrayillumination)
if xrayillumination
    if prompt
        fprintf('Enhancing contrast...\n');
    end

    lower_bound = parameters.lowerbound;
    upper_bound = parameters.higherbound;

    if lower_bound >= upper_bound
        error('lower_bound must be less than upper_bound.');
    end

    adjustedFrames = zeros(size(preprocessedFrames), 'uint8');
    scale = 255 / (upper_bound - lower_bound);

    for i = 1:numFrames
        frame = preprocessedFrames(:,:,i);
        adjustedFrame = zeros(size(frame), 'uint8');

        adjustedFrame(frame <  lower_bound) = 0;
        adjustedFrame(frame >  upper_bound) = 255;

        mask = (frame >= lower_bound) & (frame <= upper_bound);
        adjustedFrame(mask) = uint8((double(frame(mask)) - lower_bound) * scale);

        adjustedFrames(:,:,i) = adjustedFrame;
    end

    preprocessedFrames = adjustedFrames;
end

%% Step 6: Save processed video (optional)
if savevideopreprocess
    if prompt
        fprintf('Saving the processed video...\n');
    end

    outputVideoPath = fullfile(folderSaveDir, [filename '_preprocessed.avi']);
    writerObj = VideoWriter(outputVideoPath);
    open(writerObj);

    for i = 1:numFrames
        writeVideo(writerObj, preprocessedFrames(:,:,i));
    end

    close(writerObj);

    if prompt
        fprintf('Processed video saved to %s\n', outputVideoPath);
    end
end

end

    %% Nested Function: Update Frame Callback
    function updateFrame(frameNumber, hImage, frames, hText)
        
        % Callback function to update the displayed frame based on slider value
        set(hImage, 'CData', frames(:,:,frameNumber));
        set(hText, 'String', sprintf('Frame %d', frameNumber));
    end
