function [Q,M] = projectBC(Q,BC,G,H,opts)

%% PROJECTBC.m
%
% M = PROJECTBC(BC,Nleg,Mleg,DERORD)
%
% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/10/2015
% Last Modified:    06/04/2016
% ----------------------------------------------------------------------- %

%% CODE

M = bcProjMat(BC,G,H,opts);
M = spblkdiag(1,M);
Q = M.'*(Q*M);
Q = (Q+Q.')./2;             % symmetrize again!


end


% ---------------------------------------------------------------------------- %
%                           NESTED FUNCTION
% ---------------------------------------------------------------------------- %
function M = bcProjMat(BC,G,H,opts)

% Construct matrix to project onto space of functions that satisty the
% boundary conditions specified by the matrices in the cell BC. Assume one
% of the entries of BC is nonzero. Use full matrices with the "null"
% command to avoid errors when orthogonal basis is required.


            
    % Extract data
    A = BC{1}*G;
    B = BC{2}*H;

    % Set type of null space projection
    if strcmpi(opts.BCprojectorBasis,'orth')
        projtype=[];
        cleanTol = 0;
    else
        projtype='r';
        cleanTol = 1e-15;
    end

    % Compute projector
    if opts.rigorous
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use projection for rigorous case

        nA = size(A,2);
        nB = size(B,2);
        
        if isZero(B) && ~isZero(A)
            P1 = null(full(A),projtype);
            M = spblkdiag(P1,speye(nB));

        elseif isZero(A) && ~isZero(B)
            P2 = null(full(B),projtype);
            M = spblkdiag(speye(nA),P2);

        elseif ~isZero(A) && ~isZero(B)
            z = ~any(A,2);      % logical indices of zero rows in A
            Abar = A(~z,:);
            B0 = B(z,:);
            Bbar = B(~z,:);

            % First null space projection
            if ~isempty(z)
                P2 = null(full(B0),projtype);
            else
                P2 = speye(nB);
            end

            % Second projection
            if ~isempty(P2) && ~isZero(Bbar)
                % First projection had non-trivial null space
                P1 = null( null(full(Bbar*P2),projtype)'*Abar, projtype);
            elseif ~isempty(P2)
                P1 = null(full(Abar), projtype);
            else
                % BC is simply Abar*variables = 0, and set P2=0
                P1 = null(full(Abar), projtype);
                P2 = zeros(nB);
            end

            % Set output
            M = spblkdiag(sparse(P1),sparse(P2));

        else
            M=1;

        end

    elseif ~opts.rigorous
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simpler projection if not rigorous!
        % Does not care about zero columns...
        M = A + B;
        M = sparse(null(full(M), projtype));

    else
        M = sparse(null(full([A,B]), projtype));
        
    end

    % Clean
%     M(abs(M)<=cleanTol) = 0;
    
    
% END FUNCTION
end