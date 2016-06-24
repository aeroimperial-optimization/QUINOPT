function H = expandHighOrdBndVals(Nleg,Mleg,DERORD,MAXDER)

%% EXPANDHIGHORDBNDVALS.m Expand boundary values of high order derivatives 
% (outer relaxation only!)
%
% H = EXPANDHIGHORDBNDVALS(Nleg,Mleg,DERORD,MAXDER) expands the boundary values
% of the derivatives of order from DERORD to MAXDER using the Legendre
% coefficients of the derivative of order DERORD.
% *NOTE*: this is only valid when we are dealing with an OUTER APPROXIMATION,
% meaning that the functions in the integrand inequality are POLYNOMIALS of
% degree Nleg.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    21/06/2016
% Last Modified:    21/06/2016
% ----------------------------------------------------------------------- %


ndvars = length(DERORD);    % number of dependent variables
H = [];                     % initialise empty matrix
for i = 1:ndvars
    Hblk = findH(Nleg,Mleg,DERORD(i),MAXDER(i));  
    H = spblkdiag(H,Hblk);
end

%% Nested Funtion - find H
% ----------------------------------------------------------------------- %
% FIND MATRIX H TO REPRESENT BOUNDARY VALUES OF A SINGLE FUNCTION
% ----------------------------------------------------------------------- %

    function Hblk = findH(Nleg,Mleg,k,l)
        
        % Initialize
        Hblk = zeros(2*(l-k+1),Nleg-k+1);
        D = speye(Nleg-k+1);
        
        % Derivative of order k
        Hblk(1,:) = (-1).^(0:Nleg-k);
        Hblk(2,:) = 1;
        
        % Loop over derivatives of order k-1,...,l
        for j = 1:l-k
            Np1 = (Nleg-k-j) +1;     % degree of derivative of order k-j, plus 1
            v = 1:2:2*Np1-1;
            M = spdiags(ones(Np1,ceil(Np1/2)),0:2:Np1-1,Np1,Np1);
            M = spdiags(v(:),0,Np1,Np1)*[sparse(Np1,1), M];
            D = M*D;
            E = [(-1).^(0:Np1-1); ones(1,Np1)];
            Hblk([2*j+1,2*j+2],:) = E*D;  
        end
        
        % Pad with zeros
        Hblk = [sparse(2*(l-k+1),k), Hblk, sparse(2*(l-k+1),Mleg+2*k)];

    end
% ----------------------------------------------------------------------- %
end