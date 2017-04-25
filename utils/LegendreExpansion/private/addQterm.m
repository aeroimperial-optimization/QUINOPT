function [T,S,LMI,auxVars] = addQterm(T,S,LMI,auxVars,Nleg,Mleg,degp,pcoef,nnzIdx,dvar,ALPHA,BETA,DERORD)

%% addQterm.m
%
% [T,LMIs] = addQterm(T,LMIs,Nleg,Mleg,P,dvar,ALPHA,BETA,DERORD) computes
%   the relaxation of the Q term of the Legendre expansion.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/10/2015
% Last Modified:    24/02/2016
% ----------------------------------------------------------------------- %

%% CODE

Ka = DERORD(dvar(1));   % maximum derivative of first variable
Kb = DERORD(dvar(2));   % maximum derivative of second variable

% Find the appropriate portion of the matrix T to add terms in
TpartU = cumsum( Nleg+Mleg+1+2*DERORD); % partition indices of T (upper limit)
TpartL = [1, TpartU(1:end-1)+1];        % partition indices of T (lower limit)
row = TpartL(dvar(1)):TpartU(dvar(1));
col = TpartL(dvar(2)):TpartU(dvar(2));


if ALPHA==Ka && BETA==Kb && degp > 0 
    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
    % BOTH HIGHEST ORDER DERIVATIVES
    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
    % Need to add auxiliary variables
    
    if dvar(1)==dvar(2)
        
        Nk = Nleg+Mleg+Ka;
        
        % Auxiliary variables
        R = sdpvar(degp);
        Sigma = sdpvar(1);
        
        % Set integral of Leg polys
        Y = legendreTripleProduct(nnzIdx-1,Nk+1-degp,Nk,Nk+1,Nk+degp);
        H = pcoef(1).*Y{1};
        for j = 2:length(pcoef)
            H = H + pcoef(j).*Y{j};
        end
        
        % Auxiliary LMI
        % <------
        % Old code - not rescaled
        % LMI{end+1,1} = [R, H; H', Sigma.*spdiags(2./(2*(Nk+1:Nk+degp)'+1),0,degp,degp)];
        % ------>
        LMI{end+1,1} = [R, H; H', Sigma.*speye(degp)];
        auxVars = [auxVars; {R; Sigma}];
        
        % Set
        if isa(S,'sdpvar') || isa(Sigma,'sdpvar')
            S = sdpvarAddInPlace(S,-Sigma,dvar(1),dvar(2));
        else
            S(dvar(1),dvar(2)) = S(dvar(1),dvar(2)) - Sigma;
        end
        
        ind = row(Ka+Nk+2-degp:end);
        if isa(T,'sdpvar') || isa(R,'sdpvar')
            T = sdpvarAddInPlace(T,-R,ind,ind);
        else
            T(ind,ind) = T(ind,ind) - R;
        end
        return
        
        
    elseif Ka-Kb+degp <= 0
        % Inf dim term with Vk vanishes
        
        nMin = Nleg+Mleg+Ka+1-degp;
        nMax = Nleg+Mleg+Kb;
        mMin = Nleg+Mleg+Ka+1;
        mMax = Nleg+Mleg+Kb+degp;
        
        % Auxiliary variables
        R = sdpvar(Kb-Ka+degp);
        Sigma = sdpvar(1);
        
        % Set integral of Leg polys
        Y = legendreTripleProduct(nnzIdx-1,nMin,nMax,mMin,mMax);
        H = pcoef(1).*(0.5.*Y{1});
        for j = 2:length(pcoef)
            H = H + pcoef(j).*(0.5.*Y{j});
        end
        
        % Auxiliary LMI
        % <------
        % Old code - not rescaled
        % LMI{end+1,1} = [R, H; H', Sigma.*spdiags(2./(2*(mMin:mMax)'+1),0,Kb-Ka+degp,Kb-Ka+degp)];
        % ------>
        LMI{end+1,1} = [R, H; H', Sigma.*speye(Kb-Ka+degp)];
        auxVars = [auxVars; {R; Sigma}];
        
        % Set (use indices of dvar(2) for T)
        if isa(S,'sdpvar') || isa(Sigma,'sdpvar')
            S = sdpvarAddInPlace(S,-Sigma,dvar(1),dvar(1));
        else
            S(dvar(1),dvar(1)) = S(dvar(1),dvar(1)) - Sigma;
        end
        
        ind = col(Kb+Nleg+Mleg+Ka+2-degp:end);
        if isa(T,'sdpvar') || isa(R,'sdpvar')
            T = sdpvarAddInPlace(T,-R,ind,ind);
        else
            T(ind,ind) = T(ind,ind) - R;
        end
        return
        
        
    elseif Kb-Ka+degp <= 0
         % Inf dim term with Uk vanishes
         % (Same as last case but swap Ka and Kb, dvar(1) and dvar(2))
        
        nMin = Nleg+Mleg+Kb+1-degp;     % given assumption on Nleg and Mleg, always >=0!
        nMax = Nleg+Mleg+Ka;
        mMin = Nleg+Mleg+Kb+1;
        mMax = Nleg+Mleg+Ka+degp;
        
        % Auxiliary variables
        R = sdpvar(Ka-Kb+degp);
        Sigma = sdpvar(1);
        
        % Set integral of Leg polys
        Y = legendreTripleProduct(nnzIdx-1,nMin,nMax,mMin,mMax);
        H = pcoef(1).*(0.5.*Y{1});
        for j = 2:length(pcoef)
            H = H + pcoef(j).*(0.5.*Y{j});
        end
        
        % Auxiliary LMI
        % <------
        % Old code - not rescaled
        % LMI{end+1,1} = [R, H; H', Sigma.*spdiags(2./(2*(mMin:mMax)'+1),0,Kb-Ka+degp,Kb-Ka+degp)];
        % ------>
        LMI{end+1,1} = [R, H; H', Sigma.*speye(Kb-Ka+degp)];
        auxVars = [auxVars; {R; Sigma}];
        
        
        % Set (use indices of dvar(1) for T)
        if isa(S,'sdpvar') || isa(Sigma,'sdpvar')
            S = sdpvarAddInPlace(S,-Sigma,dvar(2),dvar(2));
        else
            S(dvar(2),dvar(2)) = S(dvar(2),dvar(2)) - Sigma;
        end
        
        ind = row(Kb+Nleg+Mleg+Ka+2-degp:end);
        if isa(T,'sdpvar') || isa(R,'sdpvar')
            T = sdpvarAddInPlace(T,-R,ind,ind);
        else
            T(ind,ind) = T(ind,ind) - R;
        end
        return
        
        
    else
        % Auxiliary variables
        R = sdpvar(2*degp);
        Sigma = diag(sdpvar(2,1));
        dimuu = Kb-Ka+degp; 
        dimvv = Ka-Kb+degp;
        NKa = Nleg+Mleg+Ka;
        NKb = Nleg+Mleg+Kb;
        
        % Set integral of Leg polys
        HighDer = [Ka, Kb; Kb, Ka];
        Q = cell(2,1);
        for i=1:2
            Ka = HighDer(i,1);
            Kb = HighDer(i,2);
            nMin = Nleg+Mleg+Kb+1-degp;
            nMax = Nleg+Mleg+Ka;
            mMin = Nleg+Mleg+Kb+1;
            mMax = Nleg+Mleg+Ka+degp;
            Y = legendreTripleProduct(nnzIdx-1,nMin,nMax,mMin,mMax);
            H = pcoef(1).*Y{1};
            for j = 2:length(pcoef)
                H = H + pcoef(j).*Y{j};
            end
            Q{i} = 0.5.*H;
        end
        Y = [sparse(dimvv,dimuu), Q{1}; Q{2}, sparse(dimuu,dimvv)];
        
        % Auxiliary LMI
        % <------
        % Old code - not rescaled
        % P = blkdiag(Sigma(1,1).*spdiags(2./(2*(NKa+1:NKa+degp)'+1),0,dimuu,dimuu), ...
        %             Sigma(2,2).*spdiags(2./(2*(NKb+1:NKb+degp)'+1),0,dimvv,dimvv) );
        % ------>
        P = blkdiag(Sigma(1,1).*speye(dimuu), Sigma(2,2).*speye(dimvv));
        LMI{end+1,1} = [R, Y; Y', P];
        auxVars = [auxVars; {R; Sigma}];
        
        % Set outputs
        if isa(S,'sdpvar') || isa(Sigma,'sdpvar')
            S = sdpvarAddInPlace(S,-Sigma,dvar,dvar);
        else
            S(dvar,dvar) = S(dvar,dvar) - Sigma;
        end
        
        ind1 = row(Kb+Nleg+Mleg+Ka+2-degp:end); % indices of u (first dvar)
        ind2 = col(Ka+Nleg+Mleg+Kb+2-degp:end); % indices of v (second dvar)
        ind = [ind1,ind2];
        if isa(T,'sdpvar') || isa(R,'sdpvar')
            T = sdpvarAddInPlace(T,-R,ind,ind);
        else
            T(ind,ind) = T(ind,ind) - R;
        end
        return
    end
    
    
else
    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
    % AT MOST ONE HIGHEST ORDER DERIVATIVE
    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
    % Exact representation!
    
    DerPairs = [ALPHA, BETA; BETA, ALPHA];
    HighDer = [Ka, Kb; Kb, Ka];
    Q = cell(2,1);
    
    for i=1:2
        
        ALPHA= DerPairs(i,1);   % derivative order of the finite dimensional term
        BETA = DerPairs(i,2);   % derivative order of the infinite dimensional term
        Ka = HighDer(i,1);
        Kb = HighDer(i,2);
        
        % Limits for integration matrices
        % NOTE: since have chosen Nleg >= Lp+max(Ka,Kb)-1, always have 
        % positive nMin (in fact, nMin >= max(Ka,Kb)+BETA )
        nMin = Nleg+BETA+1-degp;
        nMax = Nleg+ALPHA;
        mMin = Nleg+BETA+1;
        mMax = Nleg+ALPHA+degp;
        
        if ( nMax>=nMin ) && ( mMax>=mMin ) 
            
            % if all limits are sorted, then must compute            
            [Da,Ba] = legendreDiff(Nleg,Mleg,ALPHA,Ka,[nMin,nMax]);
            [Db,Bb] = legendreDiff(Nleg,Mleg,BETA,Kb,[mMin, mMax]);
            
            % Matrix of integral of triple products of Legendre polynomials
            % Values for l are given by nnzIdx-1
            Y = legendreTripleProduct(nnzIdx-1,nMin,nMax,mMin,mMax);
            H = pcoef(1).*([Ba';Da']*Y{1}*[Bb, Db]);
            for j = 2:length(pcoef)
                H = H + pcoef(j).*([Ba';Da']*Y{j}*[Bb, Db]);
            end
            
            % Set matrix Q
            Q{i} = H;
            
        else
            % Nothing to do - term vanishes thanks to orthogonality of
            % Legendre polyomials
            if i==1
                 Q{i} = sparse(length(row),length(col));
            elseif i==2
                 Q{i} = sparse(length(col),length(row));
            end
        end
        
    end
    
    % Set outputs
    M = Q{1} + Q{2}';
    if isa(T,'sdpvar') || isa(M,'sdpvar')
        T = sdpvarAddInPlace(T,M,row,col);
    else
        T(row,col) = T(row,col) + M;
    end
    
    
end


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
% END OF FUNCTION
end
%% END SCRIPT
