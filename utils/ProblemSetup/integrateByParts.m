function [INEQ,equalities,FLAG] = integrateByParts(INEQ)

% INTEGRATEBYPARTS.m Integrate integral inequality by parts
%
% Integrate by parts the blocks of the matrix Fi, where the size of each
% block is inferred from the list of derivative orders MAXDER. The output
% FLAG indicates if the problem is infeasible or not.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    12/04/2015
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% CODE

x = INEQ.IVAR;          % extract the independent variable for convenience
MAXDER = INEQ.MAXDER;   % extract MAXDER for convenience
FLAG = 0;               % the default for no problem

% ----------------------------------------------------------------------- %
% NUMBER OF BLOCKS AND PARTITION
% ----------------------------------------------------------------------- %
nBlk = length(MAXDER);
UFi = cumsum(MAXDER+1);     % upper indices of diagonal blocks in Fi
LFi = [1, UFi(1:end-1)+1];  % lower indices of diagonal blocks in Fi
UFb = cumsum(2*MAXDER+2);   % upper indices of diagonal blocks in Fb
LFb = [1, UFb(1:end-1)+1];  % lower indices of diagonal blocks in Fb

% ----------------------------------------------------------------------- %
% CHECK ON Fb, Fm AND Lb
% ----------------------------------------------------------------------- %
n = size(INEQ.F.Fi,1);

% Check on Fb
if isempty(INEQ.F.Fb)
    % Reinitialise to sdpvar with dummy entry
    INEQ.F = rmfield(INEQ.F,'Fb');
    INEQ.F.Fb(2*n,2*n) = x;
elseif isnumeric(INEQ.F.Fb)
    % Add to an sdpvar with dummy entry to obtain sdpvar
    A(2*n,2*n) = x;
    INEQ.F.Fb = A+INEQ.F.Fb;
end

% Check on Fm
if ~isempty(INEQ.F.Fm) && ~isZero(INEQ.F.Fm);
    intFm = 1;
else
    % Set Fm to zero sparse matrix
    intFm = 0;
    INEQ.F.Fm = spalloc(size(INEQ.F.Fb,1),n,0);
end

% Check on Li
if ~isempty(INEQ.L.Li) && ~isZero(INEQ.L.Li);
    intLi = 1;
else
    % Set Li to zero sparse matrix
    intLi = 0;
    INEQ.L.Li = spalloc(n,1,0);
end

% Check on Lb
if isempty(INEQ.L.Lb)
    % Reinitialise to sdpvar with dummy entry
    INEQ.L = rmfield(INEQ.L,'Lb');
    INEQ.L.Lb(2*n,1) = x;
elseif isnumeric(INEQ.L.Lb)
    % Add to an sdpvar with dummy entry to obtain sdpvar
    A(2*n,1) = x;
    INEQ.L.Lb = A+INEQ.L.Lb;
end

% ----------------------------------------------------------------------- %
% INTEGRATE BY PARTS DIAGONAL BLOCKS
% ----------------------------------------------------------------------- %
% Integrate by parts diagonal blocks to make them diagonal.
% Keep a list of which blocks have rows/cols to remove, and the first
% row/col to be removed. e.g. rmBlk = [2], rmfrom = [3] with block 2 having
% size 4 means that the rows/cols 3:4 of block 2 must be removed. This is
% useful information to integrate by parts the off-diagonal blocks.

keepFi = [];        % list of rows/cols of Fi to keep
rmBlk = [];         % list of blocks with row/col to remove
rmfrom = MAXDER+2;  % list of first row/col of block of Fi to remove 
                    % (this initialization means "remove none" since block 
                    % size is MAXDER+1)

for n = 1:nBlk
    
    I = LFi(n):UFi(n);          % Indices for submatrix of Fi
    J = LFb(n):UFb(n);          % Indices for submatrix of Fb
    
    % Integrate by parts if not already diagonal
    if ~isDiagonal(INEQ.F.Fi(I,I))
        [INEQ.F.Fi(I,I),INEQ.F.Fb(J,J)] = integrateByPartsFiDiagBlock(INEQ.F.Fi(I,I),INEQ.F.Fb(J,J),x);
    end 
    
    % Find last row/col to keep
    DiagEntries = diag(INEQ.F.Fi(I,I));
    if isnumeric(DiagEntries)
        lastnnz = find(DiagEntries,1,'last');
    else % sdpvar or legpoly
        lastnnz = find(any(DiagEntries),1,'last');
    end
    
    % Update list of rows/cols of Fi to keep
    keepFi = [keepFi, I(1):I(lastnnz)];
   
    if lastnnz~=numel(I)
        % 1) Update list of blocks with rows/cols to remove
        rmBlk = [rmBlk,n];
        rmfrom(n) = lastnnz+1;
        
        % 2) Integrate by parts mixed terms (if any)
        if intFm
            [INEQ.F.Fm(:,I),INEQ.F.Fb(:,J)] = integrateByPartsFm(INEQ.F.Fm(:,I),INEQ.F.Fb(:,J),x,rmfrom(n));
        end
    end
   
end


% ----------------------------------------------------------------------- %
% INTEGRATE BY PARTS OFF-DIAGONAL BLOCKS OF Fi IF NEEDED
% ----------------------------------------------------------------------- %
equalities = [];    % initialize empty list of equality constraints

if ~isempty(rmBlk) && nBlk > 1
   
    % Find block indexing (upper triangular)
    [I,J] = meshgrid(1:nBlk,1:nBlk);
    ind = any(ismember([vec(tril(I,-1)),vec(tril(J,-1))],rmBlk),2);
    I = I(ind);
    J = J(ind);
    
    % For each block
    for n = 1:length(I)
        
        % Integrate by parts 
        ri = LFi(I(n)):UFi(I(n));    % row indices for block of Fi
        ci = LFi(J(n)):UFi(J(n));    % col indices for submatrix of Fi
        rb = LFb(I(n)):UFb(I(n));    % row indices for block of Fi
        cb = LFb(J(n)):UFb(J(n));    % col indices for submatrix of Fi
        
        [INEQ.F.Fi(ri,ci),INEQ.F.Fb(rb,cb),neweq,FLAG] = ...
            integrateByPartsFiOffDiagBlock(INEQ.F.Fi(ri,ci),INEQ.F.Fb(rb,cb),x,rmfrom(I(n)),rmfrom(J(n)));
        
        % Update list if all good
        if FLAG==1
            return
        elseif ~isempty(neweq)
            equalities = [equalities; neweq];
        end
        
    end
end



% ----------------------------------------------------------------------- %
% INTEGRATE BY PARTS LINEAR TERMS
% ----------------------------------------------------------------------- %
for n = 1:nBlk
    
    I = LFi(n):UFi(n);          % Indices for submatrix of Li
    J = LFb(n):UFb(n);          % Indices for submatrix of Lb
    if intLi
        [INEQ.L.Li(I),INEQ.L.Lb(J)] = integrateByPartsLi(INEQ.L.Li(I),INEQ.L.Lb(J),x);
    end
    
end
% ----------------------------------------------------------------------- %
% SET OUTPUTS
% ----------------------------------------------------------------------- %
INEQ.DERORD = rmfrom - 2;        % new highest order derivatives for each block
INEQ.F.Fi = INEQ.F.Fi(keepFi,keepFi);
INEQ.F.Fm = INEQ.F.Fm(:,keepFi);
INEQ.F.Fb = replace(INEQ.F.Fb,x,0);       % Eliminate fake dependence on x from Fb
INEQ.L.Li = INEQ.L.Li(keepFi);
INEQ.L.Lb = replace(INEQ.L.Lb,x,0);       % Eliminate fake dependence on x from Lb

% END CODE
end
