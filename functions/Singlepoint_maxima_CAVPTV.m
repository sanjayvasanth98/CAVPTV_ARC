function Single_point_maxima = Singlepoint_maxima_CAVPTV(parameters, folderSaveDir, fileName, inputType, videoFile, preprocessedFrames)
%---------------------------------------------------------------------------%            
            % Single Point Maxima Detection Function for CAVPTV %

    % Input Arguments:
    %   parameters: struct containing necessary parameters
    %   folderSaveDir: directory to save results
    %   fileName: name of the video or image file (without extension)
    %   inputType: 'Images' or 'Video'
    %   videoFile: path to video file or folder containing images
    %   preprocessedFrames: if not empty, contains all frames as a 3D array
    %
    % Outputs:
    %   Single_point_maxima: struct containing detected centroids and number of particles
    %       detectedCentroids_SP: cell array of detected centroids per frame
    %       numParticles_SP: array of number of maxima per frame
%---------------------------------------------------------------------------%

% Extract parameters
tolerance = parameters.tolerance;              % Adjust as needed
strict = logical(parameters.strict);           % Use strict mode or not
threshold = parameters.threshold;              % Minimum value for maxima to be considered (-Inf for no threshold)
outputType = parameters.outputType;            % Output format, e.g., 'SINGLE_POINTS'
excludeOnEdges = logical(parameters.excludeOnEdges);  % Exclude maxima at image edges
displayFigures = logical(parameters.displayFigures);  % Control displaying figures
saveMaximaVideo = logical(parameters.saveMaximaVideo); % Save maxima video or not
fileExtension = parameters.fileExtension;      % E.g., '*.tif'
lightbackground = logical(parameters.lightbackground);

% Initialize outputs
Single_point_maxima = struct();

% Set up figure visibility
if displayFigures
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');
end

% Initialize variables to store results
detectedCentroids_SP = {};
numParticles_SP = [];

% Handle preprocessedFrames
if ~isempty(preprocessedFrames)
    if ndims(preprocessedFrames) == 3
        numFrames = size(preprocessedFrames, 3);
    elseif ndims(preprocessedFrames) == 4
        numFrames = size(preprocessedFrames, 4);
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
    elseif strcmpi(inputType, 'Video')
        % videoFile is the path to the video file
        videoObj = VideoReader(videoFile);
        numFrames = floor(videoObj.Duration * videoObj.FrameRate);
    else
        error('Invalid inputType. Must be ''Images'' or ''Video''.');
    end
end

% Preallocate maximaStack
if ~isempty(preprocessedFrames)
    if ndims(preprocessedFrames) == 3
        maximaStack = false(size(preprocessedFrames));
    elseif ndims(preprocessedFrames) == 4
        maximaStack = false(size(preprocessedFrames,1), size(preprocessedFrames,2), numFrames);
    end
else
    % We need to read one frame to get the size
    if strcmpi(inputType, 'Images')
        A = imread(fullfile(folderPath, fileList(1).name));
    elseif strcmpi(inputType, 'Video')
        A = read(videoObj, 1);
    end
    maximaStack = false(size(A,1), size(A,2), numFrames);
end

% Process all frames
for i = 1:numFrames
    % Read each frame
    if ~isempty(preprocessedFrames)
        if ndims(preprocessedFrames) == 3
            ip = preprocessedFrames(:,:,i);
        elseif ndims(preprocessedFrames) == 4
            ip = preprocessedFrames(:,:,:,i);
        end
    elseif strcmpi(inputType, 'Images')
        ip = imread(fullfile(folderPath, fileList(i).name));
    elseif strcmpi(inputType, 'Video')
        ip = read(videoObj, i);
    end

    % Call the findMaxima_CAVPTV function
    maximaBinary = findMaxima_CAVPTV(ip, tolerance, strict, threshold, lightbackground, outputType, excludeOnEdges);

    % Store maxima binary image
    maximaStack(:,:,i) = maximaBinary;

    % Count the number of maxima
    numMaxima = sum(maximaBinary(:));
    numParticles_SP(i) = numMaxima;
    %fprintf('Frame %d: Number of maxima: %d\n', i, numMaxima);

    % Get coordinates of maxima
    [yMaxima, xMaxima] = find(maximaBinary);
    detectedCentroids_SP{i} = [xMaxima, yMaxima];
end

% Optionally display figures
if displayFigures
    % --- Slider for Maxima Binary Stack ---
    % Create a figure with a slider to navigate through the maxima binary frames
    hFig5 = figure('Name', 'Maxima Binary Frames');
    hAx5 = axes('Parent', hFig5);
    hImage5 = imshow(maximaStack(:,:,1)); title('Maxima Binary Frames');

    % Add a slider for frame navigation
    hSlider5 = uicontrol('Style', 'slider', 'Min', 1, 'Max', numFrames, 'Value', 1, ...
        'Units', 'normalized', 'Position', [0.2 0.01 0.6 0.05], ...
        'SliderStep', [1/(numFrames-1), 10/(numFrames-1)]);

    % Add a text label to show the current frame number
    hText5 = uicontrol('Style', 'text', 'String', 'Frame 1', 'Units', 'normalized', ...
        'Position', [0.01 0.01 0.15 0.05]);

    % Update function for the slider
    set(hSlider5, 'Callback', @(src, event) updateFrame(round(get(src, 'Value')), hImage5, maximaStack, hText5));
end

% Save Maxima Binary Frames as Grayscale AVI File
if saveMaximaVideo
    % Specify the output video filename
    outputVideoFile = fullfile(folderSaveDir, [fileName, '_maxima.avi']);

    % Create a VideoWriter object with grayscale format
    outputVideo = VideoWriter(outputVideoFile, 'Grayscale AVI');

    % Set the frame rate (optional, can set it to match the input video)
    if exist('videoObj', 'var')
        outputVideo.FrameRate = videoObj.FrameRate;
    else
        outputVideo.FrameRate = 30;  % Default frame rate
    end

    % Open the VideoWriter object
    open(outputVideo);

    % Write each frame to the video file
    for i = 1:numFrames
        % Get the current binary frame from maximaStack
        frame = maximaStack(:,:,i);

        % Convert logical frame to uint8 format (0 or 255)
        frame_uint8 = uint8(frame) * 255;

        % Write the frame to the video
        writeVideo(outputVideo, frame_uint8);
    end

    % Close the VideoWriter object
    close(outputVideo);

    disp(['Maxima binary frames video saved to ', outputVideoFile]);
end

% Outputs
Single_point_maxima.detectedCentroids_SP = detectedCentroids_SP;
Single_point_maxima.numParticles_SP = numParticles_SP;

% Nested Helper Functions
function updateFrame(frameNumber, hImage, processedFrames, hText)
    % Update the image and text based on the slider value
    set(hImage, 'CData', processedFrames(:,:,frameNumber));
    set(hText, 'String', sprintf('Frame %d', frameNumber));
end

end
