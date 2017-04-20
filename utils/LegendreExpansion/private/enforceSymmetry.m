function [Q,INEQ] = enforceSymmetry(Q,INEQ)

% ENFORCESYMMETRY.m     Enforce symmetry of dependent variables
%
% Q = ENFORCESYMMETRY(Q,INEQ) enforces the symmetry of the dependent variables
%           by setting to zero the rows/cols of Q corresponding to the even/odd
%           Legendre expansion coefficients.


% Extract useful stuff
Nleg = INEQ.Nleg;
Mleg = INEQ.Mleg;
DERORD = INEQ.DERORD;
DVAR_SYMM = INEQ.DVAR_SYMM;

% Partitioning of Q, plus add 1 to take into account linear terms
QpartU = 1+cumsum( Nleg+Mleg+1+2*DERORD ); % upper limit
QpartL = [2, QpartU(1:end-1)+1];           % lower limit

% Loop over variables with symmetry
varsID = find(DVAR_SYMM>0);
for i = varsID
    
    % Find if even/odd symmetry
    K_DER_SYMM = (-1)^(DERORD(i) + DVAR_SYMM(i));
    
    
    if K_DER_SYMM==-1
        % Keep odd coefficients only
       
        % Modify Q: set to zero even coefficients
        % Take into account that first 2*DERORD(i) entries are boundary values, not
        % Legendre coeffs
        ind = QpartL(i)+2*DERORD(i):2:QpartU(i);
        Q(ind,:) = 0;
        Q(:,ind) = 0;
        
        % Add boundary conditions - needed for derivatives up to MAXDER
        INEQ.BC = addSymmetyBCs(i,DERORD,INEQ.MAXDER,INEQ.BC,DVAR_SYMM(i));
        
    elseif K_DER_SYMM==1
        % Keep even coefficients
        
        % Modify Q: set to zero odd coefficients
        % Take into account that first 2*DERORD(i) entries are boundary values, not
        % Legendre coeffs
        ind = QpartL(i)+2*DERORD(i)+1:2:QpartU(i);
        Q(ind,:) = 0;
        Q(:,ind) = 0;
        
        % Add boundary conditions - needed for derivatives up to MAXDER
        INEQ.BC = addSymmetyBCs(i,DERORD,INEQ.MAXDER,INEQ.BC,DVAR_SYMM(i));
        
    else
        error('Woops, something went wrong.')
        
    end
    
end

% Remove linearly dependent BCs
litol = 1e-10;
INEQ.BC = lirows(INEQ.BC,litol);

% END FUNCTION
end


% --------------------------------------------------------------------------- %
% Nested function
% --------------------------------------------------------------------------- %
function BC = addSymmetyBCs(DVAR,DERORD,MAXDER,BC,SYMM)
% Add boundary conditions from symmetry
    
    % Extract BC matrices that already exist
    A = BC{1};
    B = BC{2};
    [nA,mA] = size(A);
    [nB,mB] = size(B);
    [IA,JA,VA] = find(A);
    [IB,JB,VB] = find(B);
    
    % Partition indices
    PupA = cumsum(2*DERORD);
    PlwA = [1, PupA(1:end-1)+1];
    PupB = cumsum( 2*(MAXDER-DERORD+1) );
    PlwB = [1, PupB(1:end-1)+1];
    
    % Find extra conditions
    ncond = MAXDER(DVAR)+1;
    rowsA = (nA+1:nA+DERORD(DVAR)).';
    VA = [VA; ones(DERORD(DVAR),1); (-1).^( (1:DERORD(DVAR)).' + SYMM)];
    IA = [IA; rowsA; rowsA];
    JA = [JA; (PlwA(DVAR):2:PupA(DVAR)).'; (PlwA(DVAR)+1:2:PupA(DVAR)).'];
    
    ncondB = MAXDER(DVAR)-DERORD(DVAR)+1;
    rowsB = (nB+DERORD(DVAR)+1:nB+DERORD(DVAR)+ncondB).';
    VB = [VB; ones(ncondB,1); (-1).^( (DERORD(DVAR):MAXDER(DVAR)).' + SYMM+1)];
    IB = [IB; rowsB; rowsB];
    JB = [JB; (PlwB(DVAR):2:PupB(DVAR)).'; (PlwB(DVAR)+1:2:PupB(DVAR)).'];
    
    % Set output
    BC = {sparse(IA,JA,VA,nA+ncond,mA); sparse(IB,JB,VB,nB+ncond,mB)};

end

% ----------------------------------------------------------------------- %
function [BC,idx]=lirows(BC,tol)
% Extract a linearly independent set of rows of a given matrix X
%
%    [Xsub,idx]=lirows(X)
%
% inputs:
%
%  X: The given input matrix
%  tol: A rank estimation tolerance. Default=1e-10
%
% outputs:
%
% Xsub: The extracted rows of X
% idx:  The indices (into X) of the extracted rows
    
    % Build X
    A = BC{1};
    B = BC{2};
    mA = size(A,2);
    X = [A,B];
    
    % X has no non-zeros
    if ~nnz(X)
        Xsub=[]; idx=[];
        return
    end

    if nargin<2
        tol=1e-10;
    end

    [Q,R,E] = qr(full(X'),0);   % full for compatibility with old MATLAB
    if ~isvector(R)
        diagr = abs(diag(R));
    else
        diagr = R(1);
    end

    r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
    deprows = size(X,1)-r;
    idx=sort(E(1:r));
    Xsub=X(idx,:);
    
    % Re-split BCs
    BC{1} = Xsub(:,1:mA);
    BC{2} = Xsub(:,mA+1:end);

end

% ----------------------------------------------------------------------- %
