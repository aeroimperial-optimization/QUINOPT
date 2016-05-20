function x = vec(X)

%% VEC.m Vectorize matrix

% x = VEC(X) returns a vector obained from satcking the columns of X on top
% of each other. VEC(X) is equivalent to x=X(:).

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    09/05/2016
% Last Modified:    09/05/2016
% ----------------------------------------------------------------------- %

x = X(:);

% ----------------------------------------------------------------------- %
end