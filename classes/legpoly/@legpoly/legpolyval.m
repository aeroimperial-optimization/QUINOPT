function F = legpolyval(p,x)

%% LEGPOLYVAL.m Evaluate Legendre polynomial.
% 
% F = LEGPOLYVAL(p,x) evaluates the legpoly p at x. Inputs are:
%       - p : the Legendre polynomial (class legpoly) with domain [a,b]
%       - x : the points in the interval [a,b] at which p is evaluated
%
% When p is a 1-by-1 legpoly, F is a matrix or vector of the same size as
% x containing the value of p at the corresponding value of x. When p is an
% M-by-N legpoly, F is an M-by-N cell and the entry F{i,j} is a matrix
% containing the values of p(i,j) at x.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    13/01/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

%% CODE

% if ~ismatrix(x) || ~isnumeric(x)      % Compatibility with old MATLAB - ismatrix does not exist
if  ndims(x)>2 || ~isnumeric(x)
    error('Set of points for evaluation must be a real matrix.')
end

[m,n] = size(p);
[mx,nx] = size(x);
DOM = getDomain(p);
if ~any(DOM)
    % constant polynomial!
    F = coefficients(p);
    return
end

% Rescale x , vectorize, and check all points are in [-1,1]
x = (2*x(:)-DOM(2)-DOM(1))/(DOM(2)-DOM(1));   
if any(x<-1) || any(x>1)
    error('Legendre polynomial must be evaulated ad points in its domain [%0.6g,%0.6g]',DOM(1),DOM(2));
end

% Evaluate each entry in turn
F = cell(m,n);
for i = 1:numel(p)
    
    C = p(i).coef;
    N = size(C,1);                 % number of Legendre coefficients
    
    Pnm1= zeros(size(x)); % dummy
    Pn = ones(size(x));   % P0
    pval = Pn*C(1,:);        % Initial function evaluation
    
    % If have more than 1 coefficient
    for n = 0:N-2
        Pnp1 = ((2*n+1)*x.*Pn - n*Pnm1)/(n+1);  % 3 term recurrence
        pval = pval + Pnp1*C(n+2,:);                  % Add to function evaluation
        Pnm1= Pn;
        Pn = Pnp1;
    end
    
    F{i} = reshape(pval,mx,nx);
    
end

% Set output to matrix if size 1-by-1
if numel(F)==1; F = F{1}; end

end
%% END SCRIPT
