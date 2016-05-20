function D = isDiagonal(X)

% ISDIAGONAL.m Determine whether a matrix is diagonal.
%
% ISDIAGONAL(X) returns 1 if matrix X is square and diagonal, 0 otherwise.
%    ISDIAGONAL replaces MATLAB's built-in function "isdiag" for sdpvar and 
%    legpoly objects. If X is not an sdpvar or legpoly, it is more efficient
%    to use isdiag(X).
%
% See also ISDIAG

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    05/05/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

[m,n]=size(X);
if m~=n
    D = 0;
else
    D = isZero(triu(X,1)) & isZero(tril(X,-1));
end