function x = uminus(x)

%% OVERLOADED: dvarpoly/uminus

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    16/04/2015
% Last Modified:    16/04/2016
% ----------------------------------------------------------------------- %

[m,n] = size(x);
for i = 1:m
    for j = 1:n
        x(i,j).coeff = -x(i,j).coeff;
    end
end