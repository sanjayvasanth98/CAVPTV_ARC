function [videoFiles, saveDir] = file_selection_and_save_CAVPTV(prompt)
    % Function for file selection and save path setup in PTV analysis

    if prompt == 1
        % Display options to the user
        disp('Select an option:');
        disp('1. Analyze a single file');
        disp('2. Batch processing of multiple files');

        % Prompt the user for input
        choice = input('Enter your choice (1 or 2 & press enter): ');
    else
        % Non-interactive mode: Set default choices here if needed
        choice = 2;  % For example, default to batch processing
    end

    switch choice
        case 1
            % Single file analysis
            [fileName, filePath] = uigetfile({'*.avi;*.mp4', 'Video Files (*.avi, *.mp4)'}, 'Select a video file');
            if isequal(fileName, 0)
                disp('No file selected. Exiting.');
                videoFiles = '';
                saveDir = '';
                return;
            end
            videoFiles = {fullfile(filePath, fileName)};  % Store as cell array for consistency
            saveDir = uigetdir('', 'Select folder to save the results');
            
        case 2
            % Batch processing
            [fileNames, filePath] = uigetfile({'*.avi;*.mp4', 'Video Files (*.avi, *.mp4)'}, 'Select video files', 'MultiSelect', 'on');
            if isequal(fileNames, 0)
                disp('No files selected. Exiting.');
                videoFiles = '';
                saveDir = '';
                return;
            end
            
            if ischar(fileNames)
                fileNames = {fileNames};  % Ensure it's a cell array
            end
            videoFiles = fullfile(filePath, fileNames);
            saveDir = uigetdir('', 'Select folder to save the results');
    end
end
