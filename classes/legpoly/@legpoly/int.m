function Q = int(varargin)

%% INT.m    Integrate a polynomial in the Legendre basis.
%
% P = INT(p) or P = INT(p,x) integrates the polynomial p in Legendre basis 
%       (class <a href="matlab:help('legpoly')">legpoly</a>) with respect to its independent variable x. Indefinite
%       integration is performed such that P(0)=0. If a different behaviour is
%       required, please use the function <a href="matlab:help('legpolyint')">legpolyint</a>.
%
% P = INT(p,x,a,b) computes the integral of p from a to b. The integration
%       limits a and b must be contained within the domain of definition of the
%       polynomial p, as returned by calling "getDomain(p)".
%
% See also LEGPOLY, LEGPOLYINT, @LEGPOLY/JACOBIAN, @LEGPOLY/GETDOMAIN

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

% Input checks
if nargin < 1
    error('At least one input is required.'); 
end
P = varargin{1};

if nargin < 2
    varargin{2} = getivar(P); 
end
if nargin==2
    ivar = getivar(P);
    if ~( depends(varargin{2})==depends(ivar) )
        error('You can only integrate a Legendre polynomial (class legpoly) with respect to its independent variable.');
    end
end
x = varargin{2};

if nargin < 3
    a = []; 
    b = []; 
end

if nargin == 3
    error('Two limits of integrations must be provided.'); 
end

if nargin == 4
    a = varargin{3}; 
    b = varargin{4}; 
    DOM = getDomain(P);
    if min([a,b])<DOM(1) || max([a,b])>DOM(2)
        error('Legendre polynomials can only be integrated within their domain.')
    end
end

if nargin > 4; 
    error('Too many inputs'); 
end


[m,n] = size(P);

if nargin < 3
    Q = struct([]);
elseif nargin == 4
    dummy = sdpvar(1,1);
    Q(m,n) = dummy;
end

% UGLY LOOP - but better than arrayfun?
for i = 1:numel(P)
       
    Pint = legpolyint(P(i),x);
    Pint = Pint - legpolyval(Pint,0);
    
    if nargin < 3
        % Indefinite integration
        Q(i).coef = Pint.coef;
        Q(i).ivar = Pint.ivar;
        Q(i).domn = Pint.domn;
        
    elseif nargin==4
        % Definite integration
        tmp = legpolyval(Pint,[a,b]);
        Q(i) = tmp(2)-tmp(1);
        
    end
end

if isstruct(Q); Q = legpoly(Q); end

%% END FUNCTION

%% Nested function
function Z = sdpvarAddInPlace(X,Y,varargin)

    % SDPVARADDINPLACE Add entries to submatrix of SDPVAR without calling subsasn
    % 
    % Add Y to X(I) in place. I is a linear index.

    % Get variables and basis
    Y = Y(:);
    Y_basis = getbase(Y);
    Y_vars = getvariables(Y);
    Y_length = length(Y);

    % Get input indices
    if nargin==3
        % Linear indices
        I = varargin{1};
    elseif nargin==4
        % Matrix indices
        l1 = length(varargin{1});
        l2 = length(varargin{2});

        if Y_length==l1 &&  l1==l2 
            I = sub2ind(size(X),varargin{1}(:),varargin{2}(:));
        elseif Y_length==l1*l2
            [C,R] = meshgrid(varargin{2}(:),varargin{1}(:));
            I = sub2ind(size(X),R(:),C(:));
        end
    end

    % Index within limits?
    if max(I)>numel(X) || any(I<0)
        error('Index exceeds matrix dimensions.')
    end

    % Expand Y
    [nX,mX] = size(X);
    [Iy,Jy,Vy] = find(Y_basis);
    Y_basis =  sparse(I(Iy),Jy,Vy,numel(X),1+length(Y_vars));
    Y = sdpvar(nX,mX,[],Y_vars,Y_basis);

    % Set output
    Z = X+Y;