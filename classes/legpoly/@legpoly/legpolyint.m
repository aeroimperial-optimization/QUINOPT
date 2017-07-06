function P = legpolyint(p,x,BC)

%% legpolyint.m Integrate legpoly
%
% P = LEGPOLYINT(p,x) integrates the Legendre polynomial p (class legpoly)
%       defined over the compact interval [a,b] wiuth respect to the
%       independent variable x. The constant of integration is computed 
%       such that the primitive polynomial P satisfies P(a)=0. NOTE: x is
%       only used if p is a constant polynomial; otherwise, we integrate
%       with respect to the independent variable of p. Constant polynomials
%       are integrated over the domain [-1,1] by default.
%
% P = LEGPOLYINT(p,x,BC) computes the integration constant such that P(a)=BC.
%       If p is an M-by-N polynomial, BC is either a scalar or an M-by-N
%       matrix.
%
% See also legpoly/int

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    28/02/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

%% CODE

if nargin<2;
    error('Not enough inputs.')
elseif nargin==2
    BC=zeros(size(p));
elseif nargin==3
    if numel(BC)==1
        BC = BC*ones(size(p));
    elseif any(size(p)-size(BC))
        error('Size of BC does must match the size of the polynomial.')
    end
elseif nargin>3
    error('Too many inputs.')
end

% Check integration variable
if ~isa(x,'sdpvar')
    error('Independent variable x must be an sdpvar.')
end

% Preliminaries
[m,n] = size(p);
DOM = getDomain(p);
ivar = getivar(p);
idx = find(any(ivar),1,'first');
if ~isempty(idx)
    ivar = ivar(idx);
    sf = (DOM(2)-DOM(1))/2;             % scale factor
else
    ivar = x;
    DOM = [-1,1];       % set default
    sf = 1;             % scale factor
end



% Integrate each element of p
P = struct([]);
for i=1:m
    for j=1:n
        
        % Find coefficients
        pcoef = p(i,j).coef;
        N = length(pcoef);
        v1 = 1./(2*(0:N-1)'+1);
        v2 = [1;spalloc(N-1,1,0)];
        D = spdiags([v1 v2 -v1],-1:1,N+1,N);
        Pc = sf.*(D*pcoef(:)); Pc(1) = BC(i,j) + Pc(1);
        Pc = removeTrailingZeros(Pc);
        
        % Set structure fields
        ispol = length(Pc)>1;
        P(i,j).coef = Pc;
        P(i,j).ivar = ispol*ivar;
        P(i,j).domn = ispol*DOM;
        
    end
end
P = legpoly(P);

%% END SCRIPT
end
