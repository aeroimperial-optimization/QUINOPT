function [Q,Equalities,FLAG,OriginalSize,nnzrc] = makeBlkDiag(Q,Equalities,FLAG)

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
SP = SP(nnzrc,nnzrc);     % also update sparsity pattern

% Any zero on the diagonal with nonzero off-diagonal entries?
ZD = (diag(SP)==0);
if any(ZD)
    M = tril(SP);                % only take lower triangular
    [I,J] = find(M);             % get elements to set to zero
    ind = ismember(I,find(ZD)) & (I>=J);
    ind = sub2ind(size(SP),I(ind),J(ind));
    set_to_zero = clean(Q(ind),1e-15);
    if isa(set_to_zero,'sdpvar')
        Equalities = [Equalities; set_to_zero(any(set_to_zero)>0)];
    elseif isnumeric(set_to_zero) && any(set_to_zero)
        % Woops, infeasible LMI
        FLAG = 1;
    end
    Q = Q(~ZD,~ZD);
end


% Try to detect block diagonal structure
% TO DO!

