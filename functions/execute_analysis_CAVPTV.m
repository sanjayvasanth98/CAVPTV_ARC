function [detectedCentroids_SP,Single_point_maxima,bubbleMask_all,blobMasks,bufferedcentroids]=execute_analysis_CAVPTV(sampleType, analysisType, parameters, folderSaveDir, fileName, inputType, videoFile, prompt,bgImage, maskArray, video_num)
    
    %-------------------------------------------------------------------%
    % EXECUTE_ANALYSIS Executes the appropriate analysis based on 
    % sample and analysis types.
    %
    % Syntax:
    %   execute_analysis(sampleType, analysisType, parameters, folderSaveDir, filename, inputType, videoFile)
    %
    % Inputs:
    %   sampleType    - Numeric value indicating the sample type:
    %                   1 - X-ray Samples
    %                   2 - Backlight Samples
    %                   3 - Backlight Samples Without Particles
    %   analysisType  - Numeric value indicating the analysis type:
    %                   1 - Particle Analysis
    %                   2 - Blob Analysis
    %   parameters     - Structure containing processing parameters.
    %   folderSaveDir  - Directory path to save results.
    %   filename       - Base filename for saving results (without extension).
    %   inputType      - Type of input ('Images' or 'Video').
    %   videoFile      - Path to the video file or image folder.
    %
    % Outputs:
    %   None
    %-------------------------------------------------------------------%

    switch sampleType
        case 1  % X-ray Samples
            fprintf('Processing X-ray Samples...\n');
            
            % Optional Preprocessing
            if parameters.preprocessingEnabled
                fprintf('Starting Preprocessing for X-ray Samples...\n');
                preprocessedFrames = preprocessing_CAVPTV(parameters, prompt, folderSaveDir, inputType, fileName, videoFile,bgImage);
                fprintf('Preprocessing for X-ray Samples completed.\n\n');
            else
                fprintf('Preprocessing is disabled. Skipping preprocessing step for X-ray Samples.\n\n');
            end
            
            % Particle Analysis
            if analysisType == 1
                fprintf('Starting Particle Analysis for X-ray Samples...\n');
                particles = particle_analysis_CAVPTV(parameters, folderSaveDir, fileName, inputType, videoFile, preprocessedFrames);

                fprintf('Running Single_point_maxima detection on X-ray Samples...\n');
                Single_point_maxima = Singlepoint_maxima_CAVPTV(parameters, folderSaveDir, fileName, inputType, videoFile, preprocessedFrames);

                fprintf('Buffering and identifying true particles on X-ray Samples...\n');
                Buffer_Xray = Buffer_XRay_CAVPTV(parameters, folderSaveDir, fileName, inputType, videoFile, preprocessedFrames, particles, Single_point_maxima);
                bufferedcentroids = Buffer_Xray.bufferedCentroids;

                %Blob Analysis to get blob masks and calibration values
                % [blobMasks, calibrationValues] = blob_analysis_CAVPTV(parameters, folderSaveDir, fileName, inputType, videoFile, preprocessedFrames);

                %Particle Tracking using Kalman Filter
                % trackingResults = particle_tracking_CAVPTV(bufferedCentroids, blobMasks, calibrationValues, parameters, folderSaveDir, fileName, inputType, videoFile, prompt, bgImage, maskArray);

                %Save tracking results
                % save(fullfile(folderSaveDir, [fileName, '_trackingResults.mat']), 'trackingResults');

                % Check if 'optimizedParams' are returned and update parameters
                if isfield(particles, 'optimizedParams')
                    parameters.optimizedParams = particles.optimizedParams;
                    fprintf('Optimized parameters updated from Particle Analysis.\n');
                else
                    warning('No optimizedParams returned from particle_analysis_CAVPTV.');
                end
                fprintf('Particle Analysis for X-ray Samples completed.\n\n');
            end
            
            % Blob Analysis
            if analysisType == 2
                fprintf('Starting Blob Analysis for X-ray Samples...\n');
                blob_analysis_CAVPTV(parameters, folderSaveDir, fileName, inputType, inputPath, preprocessedFrames);
                fprintf('Blob Analysis for X-ray Samples completed.\n\n');
            end
            
        case 2  % Backlight Samples
            fprintf('Processing Backlight Samples...\n');
            
            % Optional Preprocessing
            if parameters.preprocessingEnabled
                fprintf('Starting Preprocessing for Backlight Samples...\n');
                preprocessedFrames = preprocessing_CAVPTV(parameters, prompt, folderSaveDir, inputType, fileName, videoFile,bgImage);
                fprintf('Preprocessing for Backlight Samples completed.\n\n');
            else
                fprintf('Preprocessing is disabled. Skipping preprocessing step for Backlight Samples.\n\n');
            end
            
            % Single_point_maxima Analysis
            if analysisType == 1
                fprintf('Running Single_point_maxima detection on Backlight Samples...\n');
                Single_point_maxima = Singlepoint_maxima_CAVPTV(parameters, folderSaveDir, fileName, inputType, videoFile, preprocessedFrames);
                detectedCentroids_SP = Single_point_maxima.detectedCentroids_SP;
                [bubbleMask_all,blobMasks]= blob_analysismaskreturn_CAVPTV(parameters, folderSaveDir, fileName, inputType, videoFile, preprocessedFrames);
                bufferedcentroids = bufferedcentroids_backlight_CAVPTV(detectedCentroids_SP, blobMasks, parameters, folderSaveDir, fileName, maskArray);
                
                % Particle Tracking using Kalman Filter
                trackingResults = particle_tracking_CAVPTV(bufferedcentroids, bubbleMask_all, parameters, folderSaveDir, fileName, inputType, videoFile, prompt, bgImage, maskArray);

                % Save tracking results
                save(fullfile(folderSaveDir, [fileName, '_trackingResults.mat']), 'trackingResults');
                
                fprintf('Single_point_maxima detection on Backlight Samples completed.\n\n');
            end
            
            % Blob Analysis
            if analysisType == 2
                fprintf('Starting Blob Analysis for Backlight Samples...\n');
                blob_analysis_CAVPTV(parameters, folderSaveDir, fileName, inputType, inputPath, preprocessedFrames);
                fprintf('Blob Analysis for Backlight Samples completed.\n\n');
            end
            
        case 3  % Backlight Samples Without Particles
            fprintf('Processing Backlight Samples Without Particles...\n');
            
            % Optional Preprocessing
            if parameters.preprocessingEnabled
                fprintf('Starting Preprocessing for Backlight Samples Without Particles...\n');
                preprocessedFrames = preprocessing_CAVPTV(parameters, prompt, folderSaveDir, inputType, fileName, videoFile,bgImage);
                fprintf('Preprocessing for Backlight Samples Without Particles completed.\n\n');
            else
                fprintf('Preprocessing is disabled. Skipping preprocessing step for Backlight Samples Without Particles.\n\n');
            end
            
            % Only Blob Analysis is performed regardless of analysisType
            fprintf('Starting Blob Analysis for Backlight Samples Without Particles...\n');
            blob_analysis_CAVPTV(parameters, folderSaveDir, fileName, inputType, videoFile, preprocessedFrames, maskArray, video_num);
            fprintf('Blob Analysis for Backlight Samples Without Particles completed.\n\n');
            
        otherwise
            error('Unknown sample type selected. Please choose a valid sample type.');
    end
end
