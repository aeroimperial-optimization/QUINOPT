function G = expandBoundaryTerm(Nleg,Mleg,DERORD)

% EXPANDBOUNDARYTERM.m Expand boundary terms with Legendre expansions.
%
% G = EXPANDBOUNDARYTERM(Nleg,Mleg,DERORD)
%
% Construct matrix to write vector of boundary terms up to derivative Ki-1 
% in terms of the Legendre coefficients of the functions and boundary
% values at x=-1 of first Ki-1 derivatives.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    01/03/2016
% Last Modified:    05/04/2016
% ----------------------------------------------------------------------- %

ndvars = length(DERORD);    % number of dependent variables
G = [];                     % initialise empty matrix
for i = 1:ndvars
    Gblk = findG(Nleg,Mleg,DERORD(i));     % matrix G as in Lemma 2 (boundary conditions)
    G = spblkdiag(G,Gblk);
end


%% Nested Funtion - find G
% ----------------------------------------------------------------------- %
% FIND MATRIX G TO REPRESENT BOUNDARY VALUES OF A SINGLE FUNCTION
% ----------------------------------------------------------------------- %

    function Gblk = findG(Nleg,Mleg,k)
        
        % Initialise vectors for nonzero entries of G and their indices
        I = zeros(k^2+2*k,1);
        J = zeros(k^2+2*k,1);
        S = zeros(k^2+2*k,1);
        
        % Other variables: identity and zero vector
        II = speye(k);
        z = sparse(1,Nleg+Mleg+k+1);
        
        % Build entries of G
        count = 0;
        for ALPHAp1=1:k-1;
            % Derivatives from 0 to k-2
            e = II(ALPHAp1,:);
            [D,B] = legendreDiff(Nleg,Mleg,ALPHAp1,k);                  % integration matrices
            [Itmp,Jtmp,Stmp] = find([e, z; e+2.*B(1,:), 2.*D(1,:)]);    % find nonzero entries and their indices
            nInd = length(Itmp);                                        % how many nonzeros
            I(count+1:count+nInd) = Itmp + 2*(ALPHAp1-1);               % shift row indices
            J(count+1:count+nInd) = Jtmp;
            S(count+1:count+nInd) = Stmp;
            count = count+nInd;
        end
        % Derivative k-1
        [Itmp,Jtmp,Stmp] = find([II(k,:), z; II(k,:), 2, z(2:end)]);
        I(count+1:end) = Itmp + 2*k - 2;                 % shift row indices by 2*ALPHA-2
        J(count+1:end) = Jtmp;
        S(count+1:end) = Stmp;
        
        % Assemble G
        Gblk = sparse(I,J,S,2*k,Nleg+Mleg+2*k+1);
        
    end

%% END SCRIPT
end