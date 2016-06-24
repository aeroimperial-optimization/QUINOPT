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
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% Default options
options.YALMIP = [];                   % YALMIP options
options.rigorous = true;               % Rigorous relaxation (with tail term) or simple approximation
options.BCprojectorBasis = 'rref';     % set to 'orth' for orthonormal, 'rref' for rational
options.sosdeg = 6;                    % degree of polynomials for S procedure

% Assign user
if isempty(userOpts)
    return
elseif ~isstruct(userOpts)
    error('Input OPTIONS must be a structure.');
else
    fnames = fieldnames(userOpts);
    allowedNames = {'YALMIP';'rigorous';'BCprojectorBasis';'sosdeg'};
    for n=1:length(fnames)
        if any(strcmpi(fnames{n},allowedNames))
            options.(fnames{n}) = userOpts.(fnames{n});
        else
            error('Unknown field %s in input OPTIONS.',fnames{n})
        end
    end
end