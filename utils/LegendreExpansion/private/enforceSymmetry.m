function [Q,G,H,INEQ] = enforceSymmetry(Q,G,H,INEQ,opts)

% ENFORCESYMMETRY.m     Enforce symmetry of dependent variables
%
% [Q,G,H,INEQ] = ENFORCESYMMETRY(Q,G,H,INEQ,opts) enforces the symmetry of the dependent
%           variables by setting to zero the rows/cols of Q corresponding to the
%           even/odd Legendre expansion coefficients, and similarly for the
%           column of the matrix G, that represents the expansion of the
%           boundary values of the derivatives of order 0 up to INEQ.DERORD-1.


% Extract useful stuff
Nleg = INEQ.Nleg;
Mleg = INEQ.Mleg;
DERORD = INEQ.DERORD;
MAXDER = INEQ.MAXDER;
DVAR_SYMM = INEQ.DVAR_SYMM;

if opts.rigorous
    
    % Partitioning of Q, plus add 1 to take into account linear terms
    QpartU = 1+cumsum( Nleg+Mleg+1+2*DERORD ); % upper limit
    QpartL = [2, QpartU(1:end-1)+1];           % lower limit
    
    % Partitioning of G
    GpartU = cumsum( Nleg + Mleg + 1 + 2*DERORD ); % upper limit
    GpartL = [1, GpartU(1:end-1)+1];               % lower limit
    
else
    
    % Partitioning of Q, plus add 1 to take into account linear terms
    QpartU = 1+cumsum( Nleg+Mleg+1+2*(MAXDER+1) );  % upper limit
    QpartL = [2, QpartU(1:end-1)+1];                % lower limit
    
    % Partitioning of G
    GpartU = cumsum( Nleg + Mleg + 1 + 2*(MAXDER+1) ); % upper limit
    GpartL = [1, GpartU(1:end-1)+1];                   % lower limit

end

% Loop over variables with symmetry
varsID = find(DVAR_SYMM>0);
for i = varsID
    
    % Find if even/odd symmetry
    if opts.rigorous
        K_DER_SYMM = (-1)^(DERORD(i) + DVAR_SYMM(i));
    else
        K_DER_SYMM = (-1)^(MAXDER(i)+1 + DVAR_SYMM(i));
    end
    
    if K_DER_SYMM==-1
        % Keep odd coefficients only
        
        if opts.rigorous
            % Take into account that first DERORD(i) entries are boundary values, not
            % Legendre coeffs
            indQ = QpartL(i)+DERORD(i):2:QpartU(i);
            indG = GpartL(i)+DERORD(i):2:GpartU(i);

        elseif ~opts.rigorous
            % Take into account that first MAXDER(i) entries are boundary values, not
            % Legendre coeffs
            indQ = QpartL(i)+MAXDER(i)+1:2:QpartU(i);
            indG = GpartL(i)+MAXDER(i)+1:2:GpartU(i);
            H(:,indG) = 0;
        end

        Q(indQ,:) = 0;
        Q(:,indQ) = 0;
        G(:,indG) = 0;
        
        % Add boundary conditions - needed for derivatives up to MAXDER
        INEQ.BC = addSymmetyBCs(i,DERORD,MAXDER,INEQ.BC,DVAR_SYMM(i));
        
    elseif K_DER_SYMM==1
        % Keep even coefficients
        
        % Modify Q: set to zero odd coefficients

        if opts.rigorous
            % Take into account that first DERORD(i) entries are boundary values, not
            % Legendre coeffs
            indQ = QpartL(i)+DERORD(i)+1:2:QpartU(i);
            indG = GpartL(i)+DERORD(i)+1:2:GpartU(i);

        elseif ~opts.rigorous
            % Take into account that first MAXDER(i)+1 entries are boundary values, not
            % Legendre coeffs, since expanded in terms of the (l+1)-th
            % derivative
            indQ = QpartL(i)+MAXDER(i)+2:2:QpartU(i);
            indG = GpartL(i)+MAXDER(i)+2:2:GpartU(i);
            H(:,indG) = 0;
        end

        Q(indQ,:) = 0;
        Q(:,indQ) = 0;
        G(:,indG) = 0;
        

        % Add boundary conditions - needed for derivatives up to MAXDER
        INEQ.BC = addSymmetyBCs(i,DERORD,MAXDER,INEQ.BC,DVAR_SYMM(i));
        
    else
        error('Woops, something went wrong.')
        
    end
    
end

% Remove linearly dependent BCs
litol = 1e-10;
INEQ.BC = lirowsBC(INEQ.BC,litol);

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
function [BC,idx]=lirowsBC(BC,tol)
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
