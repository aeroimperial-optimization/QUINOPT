function Qm = expandMixedTerm(Nleg,Mleg,Fm,IVAR,DERORD)

% EXPANDMIXEDTERM.m
%
% Qm = EXPANDMIXEDTERM(Nleg,Mleg,Fm,IVAR,DERORD)
%
% Construct representation of mixed boundary-integral term for SDP 
% relaxation of integral inequality using Legendre expansions.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    04/03/2016
% Last Modified:    15/04/2016
% ----------------------------------------------------------------------- %

% CODE

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% SIZE OF INPUT AND PARTITION
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% dvarpart: partitioning of dependent variables in F; e.g. DERORD(i) = 4,
%           so have functions u_i, u_i', ...u_i''''. The corresponding
%           columns in Fm have indices dvarpart(i-1)+1 to dvarpart(i).

[row,col] = size(Fm);           % size of Fm
dvarpart = cumsum(DERORD+1);    % indices to partition the dependent variables
U = cumsum( Nleg+Mleg+1+2*DERORD ); % partition indices of Qm cols (upper limit)
L = [1, U(1:end-1)+1];              % partition indices of Qm cols (lower limit)

% Pre initialize Qm as dummy sdpvar (prevents bug in YALMIP)
Qm(row,U(end)) = IVAR;

for i = 1:col
    
    % For each column
    Fcol = Fm(:,i);                              % select the column
    if isZero(Fcol)
        % nothing to do!
        continue
        
    else
        dvar = find(i-dvarpart<=0,1);               % which dependent variable
        Ka = DERORD(dvar);                          % maximum derivative order
        ALPHA = Ka-dvarpart(dvar)+i;                % which derivative
        
        if ALPHA==Ka
            LIMITS = [0, Nleg+Ka];
        else
            LIMITS = [];
        end
        [Da,Ba] = legendreDiff(Nleg,Mleg,ALPHA,Ka,LIMITS); % integration matrices
        
        % Select nonzero entries only
        if isnumeric(Fcol)
            ri = find(Fcol~=0);
        else
            ri = find(any(Fcol));
        end
        Fcol = Fcol(ri);                 % the nonzero entries only
        
        % Build Fhat
        Fhat = [];
        Ntot = Nleg+ALPHA;
        for k = 1:numel(ri)
            
            % Find Legendre coefficient of entry (column vector)
            if isa(Fcol(k),'legpoly')
                pdeg = degree(Fcol(k));
                pcoef = coefficients(Fcol(k));      % an SDPVAR or numeric vector
            elseif isa(Fcol(k),'sdpvar')
                pdeg = degree(Fcol(k),IVAR);
                pcoef = legbasiscoef(Fcol(k),IVAR); % an SDPVAR or numeric vector
            elseif isnumeric(Fcol(k))
                pdeg = 0;
                pcoef = Fcol(k);                    % a scalar number
            else
                error('Unknown type for entry INEQ.F.Fm(%i,%i) (not legpoly, sdpvar or numeric).',ri(k),i)
            end
            
            % Add to Fhat
            Fhat = [Fhat; pcoef'.*2./( 2.*(0:pdeg)+1 ), spalloc(1,Ntot-pdeg,0)];               
        end
        
        % Set Qm
        Qm(ri,L(dvar):U(dvar)) = Qm(ri,L(dvar):U(dvar)) + Fhat*[Ba, Da];
        
    end    
end

% Remove spurious dependence on IVAR
Qm = replace(Qm,IVAR,0);