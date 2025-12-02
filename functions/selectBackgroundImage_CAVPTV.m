function bgImages = selectBackgroundImage_CAVPTV(inputType, videoFiles, image_extension, backgroundMode)
% SELECTBACKGROUNDIMAGE_CAVPTV
% Cluster-safe version: NO prompts, NO GUI.
%
% backgroundMode:
%   1 = Provide background stack → compute average background
%   2 = Provide single averaged background image
%   3 = No background
%
% Required for mode 1 and 2:
%   You must define the background paths in the main code.

    numVideos = numel(videoFiles);
    bgImages = cell(1, numVideos);

    fprintf('Background mode selected: %d\n', backgroundMode);

    for video_num = 1:numVideos
        
        videoFile = videoFiles{video_num};
        [~, fileName, ~] = fileparts(videoFile);

        fprintf('Processing background for video %d/%d: %s\n', ...
            video_num, numVideos, fileName);

        switch backgroundMode

            % --------------------------------------------------------
            % MODE 1: Background stack → compute average
            % --------------------------------------------------------
            case 1
                fprintf('Mode 1: Computing average background from stack for %s\n', fileName);

                % You need to define bgStackFolder in main code
                bgStackFolder = sprintf('/home/USERNAME/CAVPTV/backgrounds/%s_stack/', fileName);

                stackFiles = dir(fullfile(bgStackFolder, image_extension));

                if isempty(stackFiles)
                    warning('No background stack found for %s, using empty background.', fileName);
                    bgImages{video_num} = [];
                    continue;
                end

                sumImage = 0;
                count = 0;

                for k = 1:length(stackFiles)
                    imgPath = fullfile(bgStackFolder, stackFiles(k).name);
                    I = imread(imgPath);

                    if size(I,3) == 3
                        I = rgb2gray(I);
                    end

                    sumImage = sumImage + double(I);
                    count = count + 1;
                end

                bgImages{video_num} = uint8(sumImage / count);


            % --------------------------------------------------------
            % MODE 2: Load single averaged background image
            % --------------------------------------------------------
            case 2
                fprintf('Mode 2: Loading single averaged background for %s\n', fileName);

                % You must define this file path in your main script
                bgImagePath = sprintf('/home/USERNAME/CAVPTV/backgrounds/%s_bg.tif', fileName);

                if isfile(bgImagePath)
                    I = imread(bgImagePath);
                    if size(I,3) == 3
                        I = rgb2gray(I);
                    end
                    bgImages{video_num} = I;
                else
                    warning('Background file not found for %s, using empty background.', fileName);
                    bgImages{video_num} = [];
                end


            % --------------------------------------------------------
            % MODE 3: No background
            % --------------------------------------------------------
            case 3
                fprintf('Mode 3: No background used for %s\n', fileName);
                bgImages{video_num} = [];


            otherwise
                error('Invalid backgroundMode value (must be 1, 2, or 3).');
        end
    end

    fprintf('Background processing completed.\n');
end
