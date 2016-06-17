function prd = legMulx(c,DOMAIN)

%% LEGMULX coefficients of Legendre polynomial after multiplication by x
%
% prd = LEGMULX(c) returns the coefficients of the product x*P, where P is the
% polynomial whose Legendre coefficients are given by c and x is the independent
% variable. Allow for domains different from [-1,1].
%
% NOTE: The multiplication uses the recursion relationship for Legendre
%       polynomials in the form
%
%       xP_i(x) = ((i + 1)*P_{i + 1}(x) + i*P_{i - 1}(x))/(2i + 1)
%
% Adapted from the Python package <a href
% ="https://github.com/numpy/numpy/blob/v1.10.0/numpy/polynomial/legendre.py">numpy</a>.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    17/06/2016
% Last Modified:    17/06/2016
% ----------------------------------------------------------------------- %


%% Set default domain if not specified
if nargin < 2
    DOMAIN = [-1,1];
end

% Checks and sort
if ~isnumeric(DOMAIN) || numel(DOMAIN)~=2
    error('Input DOMAIN must be a numeric vector, DOMAIN=[a,b].')
elseif DOMAIN(1)==DOMAIN(2)
    error('Input DOMAIN must be a valid domain, DOMAIN=[a,b] with a<b.')
end
DOMAIN = sort(DOMAIN);

%% Remove trailing zeros from coefficients
c = removeTrailingZeros(c(:));
degc = length(c)-1;

%% First, assume the independent variable is in [-1,1]
% Zero polynomial first
if degc == 1 && c(1) == 0
    prd = c;
else
    m = (0:degc+1)';
    a = m./(2*m-1);
    b = (m+1)./(2*m+3);
    prd = [0; a(2:end).*c] + [b(1:end-2).*c(2:end); 0; 0];
end

%% Then account for rescaling of independent variable to DOMAIN
prd = (DOMAIN(2)-DOMAIN(1))/2.*prd + (DOMAIN(2)+DOMAIN(1))/2.*[c;0];




end

