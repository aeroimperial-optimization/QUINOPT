function [G,H] = expandBoundaryVals2(Nleg,Mleg,INEQ)

% EXPANDBOUNDARYVALS2.m Expand boundary terms with Legendre expansions.
%
% G = EXPANDBOUNDARYVALS2(Nleg,Mleg,DERORD)
%
% Construct matrix to write vector of boundary terms up to derivative Li+1 
% in terms of the Legendre coefficients of the functions and boundary
% values at x=-1 of first Li derivatives. This function is used to compute OUTER
% APPROXIMATIONS ONLY.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    01/03/2016
% Last Modified:    05/04/2016
% ----------------------------------------------------------------------- %


% G - expand boundary values up to derivative Ki-1
% H - expand boundary values of derivatives Ki to Li

% initialise empty matrices
G = [];                    
H = [];
for i = 1:INEQ.ndvars
    K = INEQ.DERORD(i);
    L = INEQ.MAXDER(i);
    
    % Find Gblk
    Gblk = findGblk(Nleg,Mleg,0,K-1,L+1,0,INEQ.DVAR_SYMM(i)); % matrix G as in Lemma 2 (boundary conditions)
    G = spblkdiag(G,Gblk);
    
    % Find Hblk
    Hblk = findGblk(Nleg,Mleg,K,L,L+1,0,INEQ.DVAR_SYMM(i));
    H = spblkdiag(H,Hblk);
end


end
%% END SCRIPT