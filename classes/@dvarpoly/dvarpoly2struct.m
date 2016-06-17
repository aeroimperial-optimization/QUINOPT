function y = dvarpoly2struct(x)

%% DVARPOLY2STRUCT.m Convert dvarpoly to structure

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    16/04/2015
% Last Modified:    16/04/2016
% ----------------------------------------------------------------------- %

[m,n] = size(x);
[y(1:m,1:n).coeff] = deal(x.coeff);
[y(1:m,1:n).ivars] = deal(x.ivars);
[y(1:m,1:n).monom] = deal(x.monom);