function D = diag(X)

%% @LEGPOLY/DIAG.m
%
% Extract diagonal from square legpoly matrix, or set vector of legpoly
% polynomials as the diagonal of a square diagonal matrix. 
% NOTE: Only the main diagonal can be accessed.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

[m,n] = size(X);

if m==1 && n==1
    D = X;
    
elseif m==1 || n==1
    % Not both, altrady checked - set diagonal
    D = X(1);
    for i=2:max(m,n); D = blkdiag(D,X(i)); end
    
elseif m==n
    % but not both one - extract diagonal
    D = X(sub2ind([m,m],(1:m)',(1:m)'));
    
else
    error('X must be a vector or a square matrix.')
    
end