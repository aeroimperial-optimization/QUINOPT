function [Q,S,slacks,MatrixInequalities,AuxVars,BCproj] = expandIntegrand(INEQ,N,opts)

%% EXPANDINTEGRAND.m Expand quadratic integral inequality using Legendre expansions.
%
% [Q,S,slacks,MatrixInequalities,AuxVars,BCproj] = EXPANDINTEGRAND(INEQ,Nleg,Mleg)
%       returns the relaxation of the integral inequality specified by the
%       input structure INEQ using the Legendre expansions with parameters
%       Nleg and Mleg. The outputs are:
%
%       - Q: the finite-dimensional matrix relaxation involving the first
%            Nleg+Mleg+INEQ.DERORD Legendre coefficients
%       - S: the matrix representing the tail term
%       - slacks: list of slack variables used in the relaxation
%       - MatrixInequalities: list of matrices for additional matrix
%                             inequalities
%       - AuxVars: list of auxiliary variables (excluding slacks) used to
%                  set up the relaxation
%       - BCproj: matrix projecting the Legendre coefficients on the
%                 boundary conditions.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    11/04/2016
% Last Modified:    15/04/2016
% ----------------------------------------------------------------------- %

%% CODE


% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Compute relaxation parameters
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
k = max(INEQ.DERORD);                % maximum order of derivative
Lp = degree(INEQ.F.Fi,INEQ.IVAR);    % maximum degree of polynomials in Fi
Lm = degree(INEQ.F.Fm,INEQ.IVAR);    % max degree in Fm
if opts.rigorous;
    Mleg = Lp+k;
else
    Mleg = 0;
end

if isempty(N)
    % Compute default value
    % Consider N=(Nleg+1) Legendre coefficients, from 0 to Nleg
    Nleg = max(Lp+k-1,Lm);
elseif ~isnumeric(N) || ~isscalar(N) || N<1 || rem(N,1)~=0
    % N was provided by the user, but is not a number or empty
    error('Input N must be an positive integer.')
else
    Nleg = max([N,Lp+k-1,Lm]);        % Consider N=(Nleg+1) Legendre coefficients, from 0 to Nleg
    if Nleg>N
        fprintf(['WARNING: You requested N = %i, but it is too small. '....
            'Using the minimum value N = %i instead.\n'],N,Nleg);
    end
end


% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Expand terms
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
ndvars = length(INEQ.DERORD);
dimint = ndvars*(Nleg+Mleg+1)+2*sum(INEQ.DERORD);   % size of terms that are expanded
dimbnd = ndvars*2+2*sum(INEQ.MAXDER-INEQ.DERORD);   % size of boundary variables not expanded

isFi = ~isempty(INEQ.F.Fi) & ~isZero(INEQ.F.Fi);
isFm = ~isempty(INEQ.F.Fm) & ~isZero(INEQ.F.Fm);
isFb = ~isempty(INEQ.F.Fb) & ~isZero(INEQ.F.Fb);
isBC = ~isempty(INEQ.BC) & ~isZero(INEQ.BC);


% Expand boundary variables with Legendre series if needed for Fb, Fm or BC
if isFb || isFm || isBC
    INEQ = permuteData(INEQ);
    G = expandBoundaryTerm(Nleg,Mleg,INEQ.DERORD,opts.rigorous);
    if opts.rigorous
        H = speye(dimbnd);
        P = spblkdiag(G,H);
    else
        % Fix expansion and set coefficients to 0
        H = fixBoundaryTerm(Nleg,Mleg,INEQ.DERORD,INEQ.MAXDER);
        P = [G;H];
    end
end


% Expand quadratic part
[Q,S,slacks,MatrixInequalities,AuxVars] = ...
    expandQuadratic(INEQ,Nleg,Mleg,dimint,dimbnd,isFi,isFm,isFb,P,opts);



% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% Project onto boundary conditions
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
if isBC
    [Q,BCproj] = projectBC(Q,INEQ.BC,G,H,opts);
else
    BCproj=[];
end

%% END CODE
end