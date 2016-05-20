function setDomain(x,a,b)

%% SETDOMAIN.m Set domain of indvar object
%
% SETDOMAIN(x,a,b) sets the domain of the independent variable x to [a,b].
%       The input x should be of class indvar.
%
% See also INDVAR, @INDVAR/GETDOMAIN

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    25/04/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% Check domain
if ~isnumeric(a)||~isnumeric(b)||numel(a)~=numel(b)||numel(a)~=numel(x)||~isreal(a+b)
    error('Inputs a and b must be real numeric scalars.')
elseif any(isinf([a,b]))
    error('The domain for the independent variable must be bounded.')
end

% Set domain
qiimodel('setindvardomain',x(:),a(:),b(:));

% END CODE
end