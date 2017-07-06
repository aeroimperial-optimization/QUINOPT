function varargout = assume(varargin)


% ASSUME.m  Add assumption on dependent variable
%
% ASSUME(U,STR) adds the assumption specified by the string STR on the dependent
%               variable U (class <a href="matlab:help('depvar')">depvar</a>). Currently, suitable assumptions are:
%
%               Symmetry
%               ========
%               >> assume(U,'even');    assumes that U is symmetric with respect
%                                       to the midpoint of the domain of U.
%               >> assume(U,'odd');     assumes that U is anty-symmetric with respect
%                                       to the midpoint of the domain of U.
%               >> assume(U,'none');    assumes that U has no symmetry. Use this
%                                       command to remove previous assumptions.
%
% ASSUME(U1,STR1,U2,STR2,...) adds the assumptions specified by the character 
%               strings STR1, STR2, ..., on the dependent variables U1, U2, ...,
%               as if set by the commands assume(U1,STR1), assume(U2,STR2), ...
%
% See also DEPVAR, INDVAR

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    20/04/2017
% Last Modified:    20/04/2017
% ----------------------------------------------------------------------- %

% Check number of outputs: do not set anything
if nargout > 0
    error('Too many outputs: please call "assume(U,STR)" with no outputs.')
end
    
% Check varargin
if rem(nargin,2)
    error('Input must be pairs of dependent variable U and property string STR.')
    
elseif nargin==2
    
    u = varargin{1};
    str = varargin{2};
    
    % Check inputs
    if ~isa(u,'depvar')
        error('Input U must be a valid depvar object.')
    elseif ~isa(str,'char')
        error('Input STR must be a string.')
    end
    
    % Extract model for this variable
    mod = qiimodel('query');
    varInd = mod.DEPVARMODEL.DVARID==u.id;
    if ~any(varInd); error('Invalid dependent variable identifier.'); end
    
    % Assume symmetry?
    if any( strcmpi(str,{'none','even','odd'}) )
        qiimodel('adddepvarsymmetry',u.id,str);
    else
        warning('Unknown options specified by STR: I will ignore this directive.')
    end
    
else
    for i = 1:nargin/2
        assume(varargin{2*i-1},varargin{2*i});
    end
    
end

% END FUNCTION
end

