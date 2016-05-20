function [M,C] = polySymmMatrix(x,deg,n)

% POLYSYMMMATRIX.m Construct matrix of polynomials with variable coefficients
%
% M = POLYSYMMMATRIX(x,deg,n) construct an n-by-n symmetric matrix M of
%       sdpvar polynomials of degree deg in the variable x. The input
%       variable x must be an sdpvar or indvar object.
%
% [M,C] = POLYSYMMMATRIX(x,deg,n) also returns the coefficient of the all
%       polynomial entries of M, listed in the vector C. 

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    05/05/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

if ~isa(x,'sdpvar')
    error('Input x should be an sdpvar object.')
elseif isa(x,'indvar')
    x = sdpvar(x);
end

C = [];
for i = 1:n
    for j = i:n
     [p,pc] = polynomial(x,deg);
     M(i,j) = p;
     M(j,i) = p;
     C = [C; pc(:)];
    end
end

