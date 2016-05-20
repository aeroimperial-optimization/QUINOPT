function x = getivar(P)

%% GETIVAR.m Find independent variable of a legpoly object

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

[m,n] = size(P);
x = [P.ivar]';
x = reshape(x,m,n);