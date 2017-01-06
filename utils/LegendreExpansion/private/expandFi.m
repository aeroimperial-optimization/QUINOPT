function [Qi,S,slk,MatrixInequalities,AuxVars] = expandFi(Qi,Nleg,Mleg,Fi,IVAR,DERORD,opts)

% EXPANDFIM.m
%
% [Qi,S,slk,MatrixInequalities,AuxVars] = EXPANDFIM(Qi,Nleg,Mleg,Fi,IVAR,DERORD,opts);
%
% Compute SDP relaxation of integral-integral terms specified by the matrix
% Fi. The entries of Fi are either scalars (class double), Legendre 
% polynomials (class legpoly), or sdpvar polynomials from YALMIP.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2016
% Last Modified:    15/04/2016
% ----------------------------------------------------------------------- %

% CODE

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% SIZE OF INPUT AND PARTITION
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% ntot: number of variables, including derivatives, specified by the user
%       through DERORD
% ndvars: number of dependent variables, not including derivatives
% row,col:  total number of dependent variables, including derivatives, as
%           extracted from Fi
% dvarpart: partitioning of dependent variables in F; e.g. DERORD(i) = 4,
%           so have functions u_i, u_i', ...u_i'''', the corresponding
%           entries in Fi have indices dvarpart(i-1)+1 to dvarpart(i).

ntot = sum(DERORD+1);
ndvars = length(DERORD);        % number of dependent variables (not including derivatives)
dvarpart = cumsum(DERORD+1);    % indices to partition the dependent variables
[row,col] = size(Fi);

% Consistency checks
if row~=col;
    error('Fi must be a square cell array.');
elseif row~=ntot
    error('The dimensions of Fi and the number of variables and derivatives is not consistent.');
end


% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% INITIALIZE VARIABLES
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
if opts.rigorous
    % Initialize matrix for tail term
    dum = sdpvar(1,1);              % a dummy variable to make sure all works
    S(ndvars,ndvars) = dum;
else
    % just set to 1 for simplicity (will not be changed)
    S = 1;
end
slk.t = []; slk.pcoef = []; % Initialize empty slacks
AuxVars = {};               % Initialize empty cell for auxiliary variables
MatrixInequalities = {};    % Initialize empty cell of auxiliary LMIs

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% LOOP OVER ENTRIES OF Fi TO DECOMPOSE THEM ONE BY ONE
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Loop over upper triangular part only as Fi is upper triangular
for i = 1:ntot
    for j = i:ntot
        
        % Extract info
        p = Fi(i,j);
        if isZero(p)
            % nothing to do
            continue
        end
        
        if isa(p,'legpoly')
            pdeg = degree(p);
            pcoef = coefficients(p);      % an SDPVAR or numeric vector
        elseif isa(p,'sdpvar')
            pdeg = degree(p,IVAR);
            pcoef = legBasisCoef(p,IVAR); % an SDPVAR or numeric vector
        elseif isnumeric(p)
            pdeg = 0;
            pcoef = p;                    % a scalar number
        else
            error('Unknown type for entry INEQ.F.Fi(%i,%i) (not legpoly, sdpvar or numeric).',i,j)
        end
        
        % Select the nonzero coefficients and their index.  
        if isa(pcoef,'sdpvar')
            % overloaded "any" in yalmip returns a matrix with entry 1(0)
            % if the corresponding entry in the input matrix is nonzero(zero)
            nnzIdx = find(any(pcoef));
        elseif isnumeric(pcoef)
            % "any" does not work as above on numeric arrays!
            nnzIdx = find(pcoef~=0);
        end
        pcoef = pcoef(nnzIdx);

        
        % Find which dependent variables we are dealing with and their
        % derivative order
        dvar = [find(i-dvarpart<=0,1), find(j-dvarpart<=0,1)];
        ALPHA = DERORD(dvar(1))-dvarpart(dvar(1))+i;
        BETA =  DERORD(dvar(2))-dvarpart(dvar(2))+j;
        
        % Compute relaxation of term \int p(x)*d^ALPHA(u)*d^BETA(v) dx
        Qi = addPterm(Qi,Nleg,Mleg,pcoef,nnzIdx,dvar,ALPHA,BETA,DERORD,opts.rigorous);
        if opts.rigorous
            % the terms Q and R only appear if series is not truncated
            [Qi,S,MatrixInequalities,AuxVars] = addQterm(Qi,S,MatrixInequalities,AuxVars,Nleg,Mleg,pdeg,pcoef,nnzIdx,dvar,ALPHA,BETA,DERORD);
            [Qi,S,slk] = addRterm(Qi,S,slk,Nleg,Mleg,pdeg,pcoef,nnzIdx,IVAR,dvar,ALPHA,BETA,DERORD);
        end
        
    end
end


% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% SET OUTPUT
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
Qi = 0.5*(Qi+Qi');
if opts.rigorous; S = replace(S,dum,0); S = 0.5*(S+S'); end


% END SCRIPT
end

