function y = sdpvar(x)

% OVERLOADED: indvar/sdpvar
%
% Remove extra features of indvar and return sdpvar instead

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    25/04/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

x = struct(x);
y = x.sdpvar;