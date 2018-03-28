function [ options ] = wise_structure_settings(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the options
    options = wise_structure_core_options();
    
    Names = fieldnames(options);
    names = lower(Names);
    
    % A finite state machine to parse name-value pairs.
    if rem(nargin,2) ~= 0
        error('Arguments must occur in name-value pairs.');
    end
    
    expectval = 0;                          % start expecting a name, not a value
    i = 1;
    while i <= nargin
        arg = varargin{i};

        if ~expectval
            if ~ischar(arg)
                error(sprintf('Expected argument %d to be a string property name.', i));
            end

            lowArg = lower(arg);

            j = strmatch(lowArg,names);
            if isempty(j)                       % if no matches
                error(sprintf('Unrecognized property name ''%s''.', arg));
            elseif length(j) > 1                % if more than one match
                % Check for any exact matches (in case any names are subsets of others)
                k = strmatch(lowArg,names,'exact');
                if (length(k) == 1)
                    j = k;
                else
                    msg = sprintf('Ambiguous property name ''%s'' ', arg);
                    msg = [msg '(' deblank(Names{j(1)})];
                    for k = j(2:length(j))'
                        msg = [msg ', ' deblank(Names{k})];
                    end
                    msg = sprintf('%s).', msg);
                    error(msg);
                end
            end
            expectval = 1;                      % we expect a value next
        else
            eval(['options.' Names{j} '= arg;']);
            expectval = 0;
        end
        i = i + 1;
    end

    if expectval
        error(sprintf('Expected value for property ''%s''.', arg));
    end



function options = wise_structure_core_options   
    % termination criteria
    options.gradient_tol = 1e-4; % gradient tolerance 
    options.delta_tol = 1e-8;
    options.iter_limit = inf;
    options.time_limit = inf;
    
    % parameters for line search
    options.beta = 0.5;
    options.alpha = 1;
    options.sigma = 1e-4;
    options.max_linesearch_iter = 50;
    options.eig_tol = 1e-6;      % tolerance to check if matrix is PD
    
    % options for printing out results
    options.verbose = 1;