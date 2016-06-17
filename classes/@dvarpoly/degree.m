function d = degree(x)

%% OVERLOADED: dvarpoly/degree

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    16/04/2015
% Last Modified:    16/04/2016
% ----------------------------------------------------------------------- %

x = dvarpoly2struct(x);
d = max(arrayfun(@(z)max(sum(z.monom,2)),x(:)));