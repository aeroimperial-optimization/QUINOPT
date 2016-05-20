function C = coefficients(P,x)

%% @LEGPOLY/COEFFICIENTS.m Coefficients of Legendre polynomial
%
% C = COEFFICIENTS(P) returns the coefficients of the ppolynomial P defined
%   in the Legendre basis (class legpoly). If P is a matrix of class legpoly,
%   C is a cell array of the same size where C(i,j) contains the coefficients
%   of P(i,j). If the degree of P is 0, however, then C is a matrix, since 
%   each polynomial entry in P has only one coefficient.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

%% CODE

if numel(P)==1
    C = P.coef;
else
    [m,n] = size(P);
    coef = {P.coef};
    [C{1:m,1:n}] = deal(coef{:});
    if degree(P)==0
        C = reshape([C{:}],m,n);
    end
end


end
%% END SCRIPT
