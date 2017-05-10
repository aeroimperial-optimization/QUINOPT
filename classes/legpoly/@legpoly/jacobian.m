function J = jacobian(P,~)

%% JACOBIAN.m   Differentiate a polynomial in the Legendre basis (class legpoly). 
%
% J = JACOBIAN(P,x) computes the derivative of each entry of the matrix P of
%   polynomials in the Legendre basis (class <a
%   href="matlab:help('legpoly')">legpoly</a>) with respect to the
%   independent variable x (class <a href="matlab:help('indvar')">indvar</a>).
%
% See also LEGPOLY, INDVAR, @LEGPOLY/INT, @LEGPOLY/LEGPOLYINT

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

[m,n] = size(P);
J = struct([]);

% UGLY LOOP - but better than arrayfun?
for i = 1:m
    for j = 1:n
        
        N = degree(P(i,j));
        
        if N == 0
            % differentiating a constant!
            J(i,j).coef = 0;
            J(i,j).ivar = 0;
            J(i,j).domn = [0,0];
            
        else
            
            % Need to compute new coefficients and rescale by a factor
            % depending on the domain
            v = 1:2:2*N-1;
            M = spdiags(ones(N,ceil(N/2)),0:2:N-1,N,N);
            DOMAIN = P(i,j).domn;
            s = 2/(DOMAIN(2)-DOMAIN(1));
            
            J(i,j).coef = s.*(v(:).*(M*P(i,j).coef(2:end)));
            J(i,j).ivar = P(i,j).ivar;
            J(i,j).domn = DOMAIN;
            
        end
    end
end

J = legpoly(J);
