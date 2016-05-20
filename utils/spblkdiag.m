function A = spblkdiag(A1,A2)
% SPBLKDIAG  Sparse block diagonal concatenation.
%
% A = SPBLKDIAG(A1,A2) Given two sparse matrices A1 and A2, this function
%   constructs their block diagonal concatenation
%
%    A = | A1 0  |
%        | 0  A2 |
%
%   where A is a sparse matrix.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    05/05/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

[m1,n1] = size(A1);
[m2,n2] = size(A2);
A = [sparse(A1), sparse(m1,n2); sparse(m2,n1), sparse(A2)];