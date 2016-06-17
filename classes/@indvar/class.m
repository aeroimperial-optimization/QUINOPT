function C = class(x)

% OVERLOADED: indvar/class
% Ugly fix to have class of indvar returned as sdpvar for compatibility
% with YALMIP functions.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    25/04/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

C = 'sdpvar';
