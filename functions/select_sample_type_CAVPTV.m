function sampleType = select_sample_type_CAVPTV(parameters, prompt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SELECT_SAMPLE_TYPE_CAVPTV Prompts the user to select the type of sample
    % or uses a predefined type based on the parameters for non-interactive mode.
    %
    % Syntax:
    %   sampleType = select_sample_type_CAVPTV(parameters, prompt)
    %
    % Inputs:
    %   parameters - Structure containing various processing parameters,
    %                including 'sampleType' for non-interactive mode.
    %   prompt     - Boolean flag (1 or 0):
    %                1 - Interactive mode: Prompt the user.
    %                0 - Non-interactive mode: Use parameters.sampleType.
    %
    % Outputs:
    %   sampleType - Numeric value indicating the selected sample type:
    %                1 - X-ray Samples
    %                2 - Backlight Samples with particles
    %                3 - Backlight Samples Without Particles
    %
    % Example:
    %   % Interactive mode
    %   sampleType = select_sample_type_CAVPTV(parameters, 1);
    %
    %   % Non-interactive mode with predefined sampleType
    %   parameters.sampleType = 3;
    %   sampleType = select_sample_type_CAVPTV(parameters, 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Check if 'parameters' is a structure
    if ~isstruct(parameters)
        error('Input "parameters" must be a structure.');
    end
    
    % Set default prompt to interactive mode if not provided
    if nargin < 2
        prompt = 1;
    end
    
    if prompt == 1
        % --- Interactive Mode: Prompt the User ---
        fprintf('\nSelect the type of sample:\n');
        fprintf('1. X-ray Samples\n');
        fprintf('2. Backlight Samples with particles\n');
        fprintf('3. Backlight Samples Without Particles\n\n');
        
        % Initialize sampleType
        sampleType = [];
        
        % Loop until a valid input is received
        while isempty(sampleType) || ~ismember(sampleType, [1, 2, 3])
            % Get user input
            sampleType = input('Enter your choice (1, 2, or 3) and press Enter: ');
            
            % Validate input
            if ~ismember(sampleType, [1, 2, 3])
                fprintf('Invalid selection. Please choose 1 for X-ray Samples, 2 for Backlight Samples with particles, or 3 for Backlight Samples Without Particles.\n\n');
                sampleType = [];
            end
        end
    else
        % --- Non-Interactive Mode: Use predefined sampleType ---
        % Check if 'sampleType' exists in parameters
        if ~isfield(parameters, 'sampleType')
            error('Parameter "sampleType" is missing in the "parameters" structure for non-interactive mode.');
        end
        
        sampleType = parameters.sampleType;
        
        % Validate 'sampleType'
        if ~ismember(sampleType, [1, 2, 3])
            error('Invalid "sampleType" in parameters. Must be 1 (X-ray), 2 (Backlight with particles), or 3 (Backlight Without Particles).');
        end
        
        fprintf('\nNon-interactive mode: Using predefined sample type.\n');
    end
    
    % --- Display the Selected Sample Type ---
    switch sampleType
        case 1
            fprintf('You selected: X-ray Samples.\n');
        case 2
            fprintf('You selected: Backlight Samples with particles.\n');
        case 3
            fprintf('You selected: Backlight Samples Without Particles.\n');
    end
end
