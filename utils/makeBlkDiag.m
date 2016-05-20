function [Q,OriginalSize,nnzrc] = makeBlkDiag(Q)

%% MAKEBLKDIAG.m Make matrix block diagonal.
%
% Q = MAKEBLKDIAG(Q) makes the matrix Q block diagonal, removing any zero
% rows/columns. The original size of the matrix and the position of nonzero
% rows/comulns are also returned as additional optional outputs.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    11/04/2016
% Last Modified:    11/04/2016
% ----------------------------------------------------------------------- %

% Find sparsity pattern & original size
if isa(Q,'sdpvar')
    SP = any(Q);
elseif isnumeric(Q)
    SP = sparse(spones(Q));
end
OriginalSize = size(SP,1);

% Remove zero rows/columns
nnzrc = find(any(SP,2));  % find indices of rows/cols which have nonzero elements
Q = Q(nnzrc,nnzrc);       % only keep non-zero rows/cols
% SP = SP(nnzrc,nnzrc);     % also update sparsity pattern

% Try to detect block diagonal structure
% TO DO!

