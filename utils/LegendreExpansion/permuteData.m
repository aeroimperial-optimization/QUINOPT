function INEQ = permuteData(INEQ)

%% PERMUTEDATA.m Permute data in INEQ model
%
% Permute the rows/columns of INEQ.F.Fb, the rows of INEQ.F.Fm and the
% columns of INEQ.BC so the matrix can be partitioned into blocks
% corresponding to the following two groups of variables:
%
% 1) Boundary values of derivatives up to DERORD-1, that can be expanded
%    with Legendre series variables.
%
% 2) Boundary values of derivatives from DERORD to MAXDER, that cannot be
%    expanded and remain variables in their own right.
%
% The boundary condition matrix INEQ.BC is also partitioned into two
% submatrices for convenience.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    11/04/2015
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% INDEXING
U2 = cumsum( 2.*(INEQ.MAXDER+1) );          % partition indices of group 2 of vars (upper limit)
U1 = U2 - 2.*(INEQ.MAXDER-INEQ.DERORD+1);   % partition indices of group 1 of vars (upper limit)
L2 = U1+1;                                  % partition indices of group 2 of vars (lower limit)
L1 = [1, U2(1:end-1)+1];                    % partition indices of group 1 of vars (lower limit)

nind1 = U1-L1+1;
nind2 = U2-L2+1;

I = zeros(1,sum(nind1));
J = zeros(1,sum(nind2));

c1=0; c2=0;
for n=1:length(INEQ.DERORD)
    I(c1+1:c1+nind1(n)) = L1(n):U1(n);
    J(c2+1:c2+nind2(n)) = L2(n):U2(n);
    c1 = c1+nind1(n);
    c2 = c2+nind2(n);
end

% Permute data
if ~isempty(INEQ.F.Fb); INEQ.F.Fb = INEQ.F.Fb([I,J],[I,J]); end
if ~isempty(INEQ.F.Fm); INEQ.F.Fm = INEQ.F.Fm([I,J],:); end
if ~isempty(INEQ.BC);   INEQ.BC = {INEQ.BC(:,I); INEQ.BC(:,J)}; end

% END FUNCTION
end