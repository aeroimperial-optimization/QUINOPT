function U = triu(X,k)

% OVERLOADED: legpoly/triu

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

if nargin < 2
    k=0;
end

[m,n] = size(X);
Y = triu(ones(m,n),k);
U = Y.*X;