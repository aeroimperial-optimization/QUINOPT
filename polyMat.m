function [P,C,V] = polyMat(x,DEG,SIZE,FLAG)

% POLYMAT.m Construct matrix of polynomials with variable coefficients
%
% P = POLYMAT(x,DEG,SIZE) with SIZE=[M N] construct an M-by-N matrix P of
%       sdpvar polynomials of degree DEG in the variable x. The input
%       variable x must be an sdpvar or indvar object. P = POLYMAT(x,DEG,N)
%       can be used as a short-hand notation for P = POLYMAT(x,DEG,[N,N]).
%
% P = POLYMAT(x,DEG,SIZE,'symm') with SIZE = [N N], or its short-hand
%       syntax P = POLYMAT(x,DEG,N,'symm') construct a symmetric N-by-N
%       matrix of polynomials. When SIZE is a vector, an error is thrown if
%       SIZE(1)~=SIZE(2).
%
% [P,C,V] = POLYMAT(...) also returns the coefficient of the polynomial
%       entries of P in the cell matrix C, and a monomial basis V, such
%       that C{i,j}=coefficients(P(i,j),x) and P(i,j)=C{i,j}'*V.
%
% See also @SDPVAR/POLYNOMIAL, YALMIP/extras/COEFFICIENTS

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    12/05/2016
% Last Modified:    12/05/2016
% ----------------------------------------------------------------------- %


% Check inputs
if nargin < 3
    error('Not enough inputs.')
elseif nargin == 3
    FLAG = 'full';
elseif nargin > 4
    error('Too many inputs.')
end


% Check x
if ~isa(x,'sdpvar')
    error('Input x should be an sdpvar object.')
elseif isa(x,'indvar')
    x = sdpvar(x);
    % Check DEG
elseif ~isnumeric(DEG) || DEG<0 || rem(DEG,1)~=0
    error('Input DEG must be a non-negative integer.')
    % Check SIZE
elseif ~isnumeric(SIZE) || numel(SIZE)>2 || isempty(SIZE)
    error('Input SIZE must be a numeric vector of one or two integers.')
    % Check FLAG
elseif numel(SIZE)==2 && SIZE(1)~=SIZE(2) && strcmpi(FLAG,'symm')
    error('You cannot construct a rectangular symmetric matrix!')
end

% Set matrix size
if numel(SIZE)==1
    SIZE = [SIZE, SIZE];
end

% Construct P
C = cell(SIZE(1),SIZE(2));
for i = 1:SIZE(1)
    % Symmetric matrix case
    if strcmpi(FLAG,'symm')
        for j = i:SIZE(2)
            [q,qc,V] = polynomial(x,DEG);
            P(i,j) = q;
            P(j,i) = q;
            C{i,j} = qc(:);
            C{j,i} = qc(:);
        end
    else
    % Full matrix
        for j = 1:SIZE(2)
            [q,qc,V] = polynomial(x,DEG);
            P(i,j) = q;
            C{i,j} = qc(:);
        end
    end
end

