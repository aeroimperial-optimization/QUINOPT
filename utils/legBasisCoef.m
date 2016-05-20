function c = legBasisCoef(p,x)

%% LEGBASISCOEF.m Coefficients in Legendre basis
%
% C = LEGBASISCOEF(P,x) computes the coefficients to express the polynomial
%       P in the Legendre basis. The  Legendre decomposition is computed so
%       that if Ln(x) is the n^th Legendre polynomial, then
%
%       P(x) = C(1)*L0(x) + C(2)*L1(x) + ... + C(n+1)*Ln(x)
%
%       Inputs are a polynomial P with numeric/sdpvar coefficients and the
%       independent variable x (an sdpvar object).

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    18/02/2015
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

%% CODE

[pc,v] = coefficients(p,x);
deg = degree(v,x,1);        % degree of x associated with each coefficient (a vector!!!)
pcoefs(deg+1)=pc;           % all coefficients, padded with zeros.
n = max(deg);               % Degree of polynomial p

%% Set coefficient of Legendre polynomials up to degree n
c = (zeros(n+1,n+1)); % polynomials go from 0 to nth degree

if n==0
    c(1,end) = 1;
elseif n==1
    c(1,end) = 1;
    c(2,end-1) = 1;
else
    c(1,end) = 1;
    c(2,end-1) = 1;
    for k = 2:n
        c(k+1,:) = ((2*k-1)/(k))*circshift(c(k,:),[0,-1])-((k-1)/(k))*c(k-1,:);
    end
end
c = flipud(c');

%% Compute coefficients of Legendre decomposition
c = c\pcoefs(:);

%% END CODE
end