function s = mpower(X,Y)

%% mpower.m
%
% c = mpower(X,Y)

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/05/2017
% Last Modified:    08/05/2017
% ----------------------------------------------------------------------- %

%% CODE

% Check inputs - only positive scalars are allowed, and square matrix
if ~isnumeric(Y) || any(Y(:)<0) || any(rem(Y(:),1)) || ~isscalar(Y(:))
    error('The power must be a non-negative integer.')
elseif size(X,1)~=size(X,2)
    error('you can only take the power of square matrices.')
end


% Possibly slow!
if isscalar(X)
    s = X.^Y;
else
    s = eye(size(X));
    for m = 1:Y
        s = s*X;
    end
end

%% END SCRIPT
end