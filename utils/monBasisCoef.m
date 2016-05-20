function MC = monBasisCoef(LC)

%% MONBASISCOEF.m Convert Legendre coefficients to monomial basis coefficients
%
% MC = MONBASISCOEF(LC) converts the Legendre coefficients LC of a polynomial
%   expressed in the Legendre basis to the coefficients MC of the same 
%   polynomial in the usual monomial basis. That is, if the Legendre
%   representation of a polynomial P is
%
%   P(x) = LC(1)*L0(x) + LC(2)*L1(x) + ... + LC(n+1)*Ln(x)
%
%   where if Ln(x) is the n^th Legendre polynomial, its monomial basis 
%   representation is
%
%   P(x) = MC(1)*L0(x) + MC(2)*L1(x) + ... + MC(n+1)*Ln(x)
%
% NOTE: Accuracy may be an issue for large degree polynomials.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    18/02/2015
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %


% Set coefficient of Legendre polynomials up to degree n

n = length(LC)-1;       % degree of the polynomial
MC = (zeros(n+1,n+1)); % polynomials go from 0 to nth degree

if n==0
    MC(1,end) = 1;
elseif n==1
    MC(1,end) = 1;
    MC(2,end-1) = 1;
else
    MC(1,end) = 1;
    MC(2,end-1) = 1;
    for k = 2:n
        MC(k+1,:) = ((2*k-1)/(k))*circshift(MC(k,:),[0,-1])-((k-1)/(k))*MC(k-1,:);
    end
end
MC = flipud(MC');
MC = MC*LC;