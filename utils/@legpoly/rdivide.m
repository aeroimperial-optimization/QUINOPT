function s = rdivide(X,Y)

% OVERLOADED: rdivide
% element-wise division

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

if isnumeric(Y) && isa(X,'legpoly')
    Y = 1./Y;
    s = X.*Y;
elseif isa(X,'double') && isa(Y,'legpoly')
    error('Only division of legpoly object by a numeric scalar is supported.')
end