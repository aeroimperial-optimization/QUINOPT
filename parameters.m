function p = parameters(varargin)

% PARAMETERS.m Create symbolic optimization variables for QUINOPT
%
% p = PARAMETERS(m,n) creates an m-by-n matrix of parameters p to be used
%       as optimization variables with QUINOPT. The parameters are
%       YALMIP variables (sdpvar objects). Multiple scalar parameters can
%       be defined with the simplified syntax
%
%       PARAMETERS a b c d
%
% NOTE: this function wraps YALMIP's <a href="matlab:help('sdpvar')">sdpvar</a> function.
%
% See also SDPVAR, INDVAR, DEPVAR, QUINOPT, SETQUADINTINEQ

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    05/05/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% NOTE: Weird wrapper to SDPVAR from YALMIP. Calling
% [varargout{:}]=sdpvar(varargout{:}) does not work as sdpvar does not
% return outputs when calling "sdpvar a b c d"; rather, it assigns the sdpvars
% directly by a call to "assignin".

if nargin == 0
    p = sdpvar(1,1);
end

if ischar(varargin{1})
    % Construct YALMIP variables with given name
    for k = 1:nargin
        if ~isvarname(varargin{k})
            error('Invalid variable name.')
        else
            assignin('caller',varargin{k},eval('sdpvar(1,1);'));
        end
    end
    return
    
elseif isnumeric(varargin{1})
    % Construct YALMIP variable with given size
    if nargin == 1
        p = sdpvar(varargin{1},1);
    elseif nargin == 2
        p = sdpvar(varargin{1},varargin{2});
        
    else
        error('Too many input')
    end
    
end