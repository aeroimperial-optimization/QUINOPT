function Qm = expandLi(Nleg,Mleg,Li,IVAR,DERORD,rigorous)

% EXPANDLI.m
%
% Qm = EXPANDLI(Nleg,Mleg,Li,IVAR,DERORD)
%
% Construct representation of linear integral term for SDP 
% relaxation of integral inequality using Legendre expansions.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    17/04/2017
% Last Modified:    17/04/2017
% ----------------------------------------------------------------------- %

% CODE

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% SIZE OF INPUT AND PARTITION
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% dvarpart: partitioning of dependent variables in F; e.g. DERORD(i) = 4,
%           so have functions u_i, u_i', ...u_i''''. The corresponding
%           entries in Li have indices dvarpart(i-1)+1 to dvarpart(i).

numLi = length(Li);             % size of Li
dvarpart = cumsum(DERORD+1);    % indices to partition the dependent variables
U = cumsum( Nleg+Mleg+1+2*DERORD ); % partition indices of Qm cols (upper limit)
L = [1, U(1:end-1)+1];              % partition indices of Qm cols (lower limit)

% Pre initialize Qm as dummy sdpvar (prevents bug in YALMIP)
% This is a row vector for convenience
Qm(1,U(end)) = IVAR;

% For each element
for i = 1:numLi

    if isZero(Li(i))
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
        
        
        % Find Legendre coefficient of entry Li(i) as a column vector
        if isa(Li(i),'legpoly')
            pdeg = degree(Li(i));
            pcoef = coefficients(Li(i));      % an SDPVAR or numeric vector
        elseif isa(Li(i),'sdpvar')
            pdeg = degree(Li(i),IVAR);
            pcoef = legBasisCoef(Li(i),IVAR); % an SDPVAR or numeric vector
        elseif isnumeric(Li(i))
            pdeg = 0;
            pcoef = Li(i);                    % a scalar number
        else
            error('Unknown type for entry INEQ.L.Li(%i) (not legpoly, sdpvar or numeric).',i)
        end
        
        % Build Lhat
        Ntot = Nleg+ALPHA;
        Lhat = [pcoef'.*sqrt(2)./sqrt( 2.*(0:pdeg)+1 ), sparse(1,Ntot-pdeg)];
        
        % Set Qm
        if ~rigorous; 
            % Set entries of unused coefficients to 0.
            Da(:,Nleg+2:end) = 0;   
        end
        Qm(L(dvar):U(dvar)) = Qm(L(dvar):U(dvar)) + Lhat*[Ba, Da];
        
    end    
end

% Remove spurious dependence on IVAR & fix for outer approximation if required
Qm = replace(Qm,IVAR,0);