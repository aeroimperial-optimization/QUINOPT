function options = setQUINOPTOptions(userOpts)

%% SETQUINOPTOPTIONS.m Set user options for QUINOPT
%
% options = SETQUADINTINEQOPTIONS(userOpts) takes the user input structure
%       options and overrides the default options.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    27/04/2015
% Last Modified:    18/05/2017
% ----------------------------------------------------------------------- %

% A useful persistent variable
persistent added_rigorous
if isempty(added_rigorous); added_rigorous = false; end

% Default options
options.YALMIP           = struct([]);  % YALMIP options
options.N                = [];          % Number of Legendre coefficients to use in the expansion
options.solve            = true;        % Solve the optimization problem? (true/false, default: true)
options.method           = 'inner';     % Inner or outer approximation ('inner' or 'outer', defauls: 'inner')
options.BCprojectorBasis = 'rref';      % set to 'orth' for orthonormal, 'rref' for rational (default: 'rref')
options.sosdeg           = 6;           % degree of polynomials for S procedure
options.psdtol           = -1e-8;       % tolerance for small negative eigenvalues

% Class of each option
opcls.YALMIP           = 'struct';
opcls.solve            = {'double';'logical'};
opcls.N                = {'double';'single';'uint8';'uint16';'uint32';'uint64';'int8';'int16';'int32';'int64'};
opcls.method           = 'char';
opcls.BCprojectorBasis = 'char';
opcls.sosdeg           = 'double';
opcls.psdtol           = 'double';

% Assign user
if isempty(userOpts)
    return
    
elseif ~isstruct(userOpts)
    error('Input OPTIONS must be a structure.');
    
else
    fnames = fieldnames(userOpts);
    allowedNames = fieldnames(options);
    for n=1:length(fnames)
        % Fix for deprecated option
        if strcmpi(fnames{n},'rigorous')
            if ~userOpts.rigorous
                options.method = 'outer';
            end
        else
            if any(strcmpi(fnames{n},allowedNames))
                % Get user option
                boo = userOpts.(fnames{n});
                boocls = class(boo);
                if any(strcmpi( boocls,opcls.(fnames{n}) ))
                    options.(fnames{n}) = boo;
                else
                    error('Wrong class for field %s in OPTIONS.',fnames{n})
                end
            else
                error('Unknown field %s in input OPTIONS.',fnames{n})
            end
        end
    end
    
    % Set options.rigorous, used in many functions (but deprecated)
    if strcmpi(options.method,'inner')
        options.rigorous = true;
    elseif strcmpi(options.method,'outer')
        options.rigorous = false;
    else
        error('Option OPTIONS.method must be either ''inner'' or ''outer''.')
    end
    
end
