function p = power(x,d)

%% OVERLOADED: dvarpoly/power

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    16/04/2015
% Last Modified:    16/04/2016
% ----------------------------------------------------------------------- %

if ~isnumeric(d) || numel(x)~=1 || any(rem(d,1)) || any(d<0)
    error('When calling x.^d x must be a scalar object and d an array of non-negative integers.')
end

% Reuse code
[m,n] = size(d);
for i = 1:numel(d)
    p(i) = x^d(i);
end
p = reshape(p,m,n);