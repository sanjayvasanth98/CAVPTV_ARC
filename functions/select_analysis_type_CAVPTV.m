function analysisType = select_analysis_type_CAVPTV(parameters, prompt)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SELECT_ANALYSIS_TYPE_CAVPTV Prompts the user to select the type of analysis
    % or uses a predefined type based on the parameters for non-interactive mode.
    
    % Syntax:
    %   analysisType = select_analysis_type_CAVPTV(parameters, prompt)
    
    % Inputs:
    %   parameters - Structure containing various processing parameters,
    %                including 'analysisType' for non-interactive mode.
    %   prompt     - Boolean flag (1 or 0):
    %                1 - Interactive mode: Prompt the user.
    %                0 - Non-interactive mode: Use parameters.analysisType.
    
    % Outputs:
    %   analysisType - Numeric value indicating the selected analysis type:
    %                  1 - Particle Analysis
    %                  2 - Blob Analysis
    
    % Example:
    %   % Interactive mode
    %   analysisType = select_analysis_type_CAVPTV(parameters, 1);
    
    %   % Non-interactive mode with predefined analysisType
    %   parameters.analysisType = 2;
    %   analysisType = select_analysis_type_CAVPTV(parameters, 0);
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
        fprintf('\nSelect the type of analysis:\n');
        fprintf('1. Particle Analysis\n');
        fprintf('2. Blob Analysis\n\n');
        
        % Initialize analysisType
        analysisType = [];
        
        % Loop until a valid input is received
        while isempty(analysisType) || ~ismember(analysisType, [1, 2])
            % Get user input
            analysisType = input('Enter your choice (1 or 2) and press Enter: ');
            
            % Validate input
            if ~ismember(analysisType, [1, 2])
                fprintf('Invalid selection. Please choose 1 for Particle Analysis or 2 for Blob Analysis.\n\n');
                analysisType = [];
            end
        end
    else
        % --- Non-Interactive Mode: Use predefined analysisType ---
        % Check if 'analysisType' exists in parameters
        if ~isfield(parameters, 'analysisType')
            error('Parameter "analysisType" is missing in the "parameters" structure for non-interactive mode.');
        end
        
        analysisType = parameters.analysisType;
        
        % Validate 'analysisType'
        if ~ismember(analysisType, [1, 2])
            error('Invalid "analysisType" in parameters. Must be 1 (Particle Analysis) or 2 (Blob Analysis).');
        end
        
        fprintf('\nNon-interactive mode: Using predefined analysis type.\n');
    end
    
    % --- Display the Selected Analysis Type ---
    switch analysisType
        case 1
            fprintf('You selected: Particle Analysis.\n');
        case 2
            fprintf('You selected: Blob Analysis.\n');
    end
end