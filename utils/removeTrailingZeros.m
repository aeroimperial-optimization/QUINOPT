function x = removeTrailingZeros(x)

% remove trailing zeros from vector. If all zeros or empty, returns a
% single zero.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

% Check input
if isempty(x)
   x = 0;
end

if~isvector(x)
    error('Input x must be a vector.')
end

% Different method for sdpvars
if isnumeric(x)
    n = find(x,1,'last');
else % for sdpvars
    n = find(any(x),1,'last');
end

if ~isempty(n)
    x = x(1:n);
else
    x = 0;
end