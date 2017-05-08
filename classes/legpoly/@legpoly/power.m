function s = power(X,Y)

%% power.m
%
% c = power(X,Y)

% Note: maintains MATLAB's behaviour that 0^0 gives 1.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/05/2017
% Last Modified:    08/05/2017
% ----------------------------------------------------------------------- %

%% CODE

% Check inputs - only positive powers are allowed
if ~isnumeric(Y) || any(Y(:)<0) || any(rem(Y(:),1))
    error('Powers must be non-negative integers.')
end

% Scalar or vector?
if ~isscalar(Y) && ~any(size(X)==size(Y))
    error('Dimensions of inputs do not match')
elseif isscalar(Y)
    Y = Y.*ones(size(X));
end

% Reuse code, possibly slow!
s = legpoly(ones(size(X)));
for n = 1:numel(X)
    % Loop multiplication
    for m = 1:Y(n)
        s(n) = s(n).*X(n);
    end
end

%% END SCRIPT
end