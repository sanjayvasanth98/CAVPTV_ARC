function maskArray = maskarray_CAVPTV(videoFiles, baseSaveDir, bgImages)
    % MASKARRAY_CAVPTV Generates masks from video files or image directories and saves them.
    %
    % maskArray = maskarray_CAVPTV(videoFiles, baseSaveDir, bgImages)
    %
    % Inputs:
    %   videoFiles    - Cell array of paths to video files or directories containing images.
    %   baseSaveDir   - String specifying the base directory to save masks.
    %   bgImages      - Cell array of grayscale background images corresponding to each video.
    %                   Use empty cells ([]) if no background image is provided for a video.
    %
    % Outputs:
    %   maskArray     - 3D array containing all generated masks.
    %
    % Example:
    %   maskArray = maskarray_CAVPTV({'video1.avi', 'video2.avi'}, 'Masks', {bgImage1, []});
    
    %------------------------ Input Validation ------------------------%
    
    % Validate inputs
    if ~iscell(videoFiles) || isempty(videoFiles)
        error('videoFiles must be a non-empty cell array.');
    end
    if ~iscell(bgImages) || length(bgImages) ~= length(videoFiles)
        error('bgImages must match number of videos.');
    end

    % Mask array (store empty or mask per video)
    maskArray = cell(1, length(videoFiles));

    fprintf('Starting mask creation for %d videos...\n', length(videoFiles));

    for video_num = 1:length(videoFiles)

        videoFile = videoFiles{video_num};
        [~, fileName, ~] = fileparts(videoFile);

        fprintf('Video %d/%d: %s\n', video_num, length(videoFiles), fileName);

        % Create folder for this video
        folderSaveDir = fullfile(baseSaveDir, fileName);
        if ~exist(folderSaveDir, 'dir')
            mkdir(folderSaveDir);
        end

        % ===============================================================
        % CASE 3 — No background given → NO MASK CREATION
        % bgImages{video_num} is empty when backgroundMode = 3
        % ===============================================================
        if isempty(bgImages{video_num})
            fprintf('  No background → No mask generated.\n');
            maskArray{video_num} = [];
            continue;   % Skip everything else
        end

        % ===============================================================
        % CASE 1 or 2 — Background provided → use it directly
        % ===============================================================
        bgImage = bgImages{video_num};

        if ndims(bgImage) ~= 2
            error('Background image must be grayscale 2D.');
        end

        avgProjection = bgImage;   % background already computed

        % ===============================================================
        % CREATE MASK (ONLY when background exists)
        % ===============================================================
        try
            mask = CreateMask_CAVPTV(avgProjection);
        catch ME
            fprintf('  Error creating mask: %s\n', ME.message);
            maskArray{video_num} = [];
            continue;
        end

        % Store mask in cell array
        maskArray{video_num} = mask;

        % Save mask per video
        maskSavePath = fullfile(folderSaveDir, [fileName '_mask.mat']);
        try
            save(maskSavePath, 'mask');
            fprintf('  Mask saved: %s\n', maskSavePath);
        catch ME
            fprintf('  Error saving mask: %s\n', ME.message);
        end
    end

    % Save the full maskArray (cell array)
    maskArraySavePath = fullfile(baseSaveDir, 'maskArray.mat');
    try
        save(maskArraySavePath, 'maskArray');
        fprintf('All masks saved in: %s\n', maskArraySavePath);
    catch ME
        fprintf('Error saving maskArray: %s\n', ME.message);
    end

    fprintf('Mask creation completed.\n');
end
