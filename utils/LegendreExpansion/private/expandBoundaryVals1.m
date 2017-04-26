function [G,H] = expandBoundaryVals1(Nleg,Mleg,INEQ)

% EXPANDBOUNDARYVALS1.m Expand boundary terms with Legendre expansions.
%
% G = EXPANDBOUNDARYVALS1(Nleg,Mleg,DERORD)
%
% Construct matrix to write vector of boundary terms up to derivative Ki-1 
% in terms of the Legendre coefficients of the functions and boundary
% values at x=-1 of first Ki-1 derivatives. This function is used to compute
% INNER APPROXIMATIONS ONLY.
 
% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    01/03/2016
% Last Modified:    05/04/2016
% ----------------------------------------------------------------------- %


% First find G
ndvars = length(INEQ.DERORD);    % number of dependent variables
G = [];                     % initialise empty matrix
for i = 1:ndvars
    k = INEQ.DERORD(i);
    Gblk = findGblk(Nleg,Mleg,0,k-1,k,1,INEQ.DVAR_SYMM(i)); % matrix G as in Lemma 2 (boundary conditions)
    G = spblkdiag(G,Gblk);
end

% Then set up H
H = speye(INEQ.dimbnd);

end
%% END SCRIPT