function P = setDomain(P,DOMAIN)

% SETDOMAIN.m Set domain of legpoly
%
% P = SETDOMAIN(P,DOMAIN) sets the domain of independent variable of a 
% legpoly object. DOMAIN must be a row vector with 2 entries, DOMAIN=[a,b].

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

if ~isnumeric(DOMAIN) || numel(DOMAIN)~=2
    error('Input DOMAIN must be a numeric vector, DOMAIN=[a,b].')
elseif DOMAIN(1)==DOMAIN(2)
    error('Input DOMAIN must be a valid domain, DOMAIN=[a,b] with a<b.')
end

DOMAIN = sort(DOMAIN);
ivar = double(any(getivar(P)));     % 1 for entries with nonzero degree
[m,n]=size(ivar);
ivar = mat2cell(ivar,ones(m,1),ones(n,1));
domn = cellfun(@(x)x*DOMAIN,ivar,'uniformoutput',0);
[P(1:m,1:n).domn] = deal(domn{:});