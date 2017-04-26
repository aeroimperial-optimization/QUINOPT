function [T,S,slk] = addRterm(T,S,slk,Nleg,Mleg,degp,pcoef,nnzIdx,IVAR,dvar,ALPHA,BETA,DERORD,DVAR_SYMM)

%% addRterm.m
%
% [T,S,slk] = addRterm(T,S,slk,Nleg,Mleg,p,dvar,ALPHA,BETA,DERORD) computes
%   the relaxation of the R term of the Legendre expansion.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/10/2015
% Last Modified:    08/10/2015
% ----------------------------------------------------------------------- %

%% CODE

Ka = DERORD(dvar(1));   % maximum derivative of first variable
Kb = DERORD(dvar(2));   % maximum derivative of second variable
SYMM_ALPHA = DVAR_SYMM(dvar(1));
SYMM_BETA  = DVAR_SYMM(dvar(2));

% Find the appropriate portion of the matrix T to add terms.
% NOTE: only need the part with actual Legendre coefficients, not the 
% boundary values
TpartU = cumsum( Nleg+Mleg+1+2*DERORD); % partition indices of T (upper limit)
TpartL = [1, TpartU(1:end-1)+1];        % partition indices of T (lower limit)
row = TpartL(dvar(1)):TpartU(dvar(1)); 
col = TpartL(dvar(2)):TpartU(dvar(2));
% row = TpartL(dvar(1))+Ka:TpartU(dvar(1)); 
% col = TpartL(dvar(2))+Kb:TpartU(dvar(2));


if ALPHA==Ka && BETA==Kb
    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
    % BOTH HIGHEST ORDER DERIVATIVES
    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
    % Construct polynomial as sdpvar object
    pcoef(nnzIdx) = pcoef;  
    pcoef = monBasisCoef(pcoef);
    mons = monolist(IVAR,degp);
    
    % Assign
    S(dvar(1),dvar(2)) = S(dvar(1),dvar(2)) + pcoef(:)'*mons;
    
    
else
    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
    % AT MOST ONE HIGHEST ORDER DERIVATIVE
    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
    
    % Find the matrices
    ee = (Nleg+1)^(BETA-ALPHA);
    Za = computeZmatrix(Nleg,Mleg,ALPHA,Ka,SYMM_ALPHA);
    lambda_a = computelambda(Nleg,Mleg,ALPHA,Ka);
    if BETA==ALPHA && Kb==Ka
        Zb = Za;
        lambda_b = lambda_a;
    else
        Zb = computeZmatrix(Nleg,Mleg,BETA,Kb,SYMM_BETA);
        lambda_b = computelambda(Nleg,Mleg,BETA,Kb);
    end
    
    % Introduce the slacks (and reference them!)
    basemat = getbase(pcoef);                   % find base matrices for pcoef (YALMIP compact representation)
    issdpvar = sum(spones(basemat(:,2:end)),2); % vector with 0 if corresponding entry of pcoef is a numeric value
    idx_sdpvar = find(issdpvar);                % indices of entries of pcoef which are sdpvars
%     idx_double = find(~issdpvar);             % numeric entries, can take absolute value
    needslk = pcoef(idx_sdpvar);                % the  coefficients that need a slack
    numcoef = pcoef(~issdpvar);                 % the numeric coefficients
    
    if ( ~isempty(slk.t) ) && ( ~isempty(needslk) )
        % returns 1 if a slack already exists and slacks are needed
        % If so, find which variables that need slacks (up to a + or -)
        % already have one
%         slkExists = belongsTo(slk.pcoef,[needslk;-needslk]);
        [slkExists,whichslk] = belongsTo(needslk,[slk.pcoef;-slk.pcoef]);
        whichslk = mod(whichslk(slkExists),length(slk.pcoef));
        whichslk(whichslk==0)=length(slk.pcoef);
    else
        slkExists = 0;
    end
    
    if any(slkExists)
        % Need to put slacks in and some already exist
        idx_needslack = idx_sdpvar(~slkExists);     % indices of variables that need new slack
        told = slk.t(whichslk);                     % load old variables
        t = sdpvar(length(idx_needslack),1);                     % new slacks
        addToList = vec(pcoef(idx_needslack));                   % variables to add to list of variables with slack
        
        % Estimate infinity norm of p
        pinf = sum([t; told; abs(numcoef(:))]);
        
    elseif ~isempty(needslk)
        % Slacks are needed but do not exist
        
        t = sdpvar(length(idx_sdpvar),1);       % all variables need a slack
        addToList = vec(pcoef(idx_sdpvar));     % all variables replaced by slacks
        
        % Estimate infinity norm of p
        pinf = sum([t; abs(numcoef(:))]);
        
    else
        % slacks are not needed - easy!
        
        t = [];                     % no slacks to add!
        addToList = [];             % no variables replaced by slacks
        pinf = sum(abs(numcoef));   % Estimate infinity norm of p
        
        
    end
    
    
    % Set output (only add nonzero terms)
    if isa(T,'sdpvar') || isa(pinf,'sdpvar')
        T = sdpvarAddInPlace(T,-(pinf*ee).*spblkdiag(zeros(Ka),Za),row,row);
        T = sdpvarAddInPlace(T,-(pinf/ee).*spblkdiag(zeros(Kb),Zb),col,col);
    else
        T(row,row) = T(row,row) - (pinf*ee).*spblkdiag( zeros(Ka),Za );
        T(col,col) = T(col,col) - (pinf/ee).*spblkdiag( zeros(Kb),Zb );
    end

    S(dvar(1),dvar(1)) = S(dvar(1),dvar(1)) - pinf*(ee*lambda_a);
    S(dvar(2),dvar(2)) = S(dvar(2),dvar(2)) - pinf*(lambda_b/ee);
    
    slk.t = [slk.t; t(:)];
    slk.pcoef = [slk.pcoef; addToList];
    
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
% END OF MAIN FUNCTION
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
end

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
% NESTED FUNCTIONS
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%

% Find Z
function  Z = computeZmatrix(Nleg,Mleg,ALPHA,K,SYMM)
    Z = computeHmatrix(Nleg,Mleg,ALPHA,K,SYMM);
    if ALPHA < K
        T = computeTmatrix(Nleg,Mleg,ALPHA+1,K,SYMM);
        Z = Z + T;
        for i=1:K-ALPHA-1
            T = computeTmatrix(Nleg,Mleg,ALPHA+1+i,K,SYMM);
            omega = computeomega(Nleg,Mleg,ALPHA+1:ALPHA+i);
            Z = Z + prod(omega).*T;
        end
    end
end

% Find lambda
function  lambda = computelambda(Nleg,Mleg,ALPHA,Ka)
    if ALPHA < Ka 
       % standard case
       omega = computeomega(Nleg,Mleg,ALPHA+1:Ka);
       lambda = 0.5*prod( omega ) ;
    elseif ALPHA==Ka
        % easy case
        lambda = 0.5;
    else
        error('ALPHA cannot be larger than the highest-order derivative.')
    end
end

% Find H
function  H = computeHmatrix(Nleg,Mleg,ALPHA,K,SYMM)
    D = legendreDiff(Nleg,Mleg,ALPHA,K,[Nleg+ALPHA+1,Nleg+Mleg+ALPHA],SYMM);
    % Old code - no rescaling
    % <-----
    % v = 1./( 2*( Nleg+ALPHA+1:Nleg+Mleg+ALPHA )' + 1 );
    % dim = length(v);
    % H = D'*spdiags(v,0,dim,dim)*D;
    % ---->
    H = (D.'*D)./2;
end

% Find T
function  T = computeTmatrix(Nleg,Mleg,eta,K,SYMM)
    if eta < K+1
        % Interesting case
        D = legendreDiff(Nleg,Mleg,eta,K,[Nleg+Mleg+eta-1,Nleg+Mleg+eta],SYMM);
        % Old code - no rescaling
        % <-----
        % N = [2/( 2*(Nleg+Mleg+eta)-1 )^2/(2*(Nleg+Mleg+eta)+1); ...
        %      2/( 2*(Nleg+Mleg+eta)+1 )^2/(2*(Nleg+Mleg+eta)+3)];
        % ---->
        N = [1/( 2*(Nleg+Mleg+eta)-1 )/(2*(Nleg+Mleg+eta)+1); ...
             1/( 2*(Nleg+Mleg+eta)+1 )/(2*(Nleg+Mleg+eta)+3)];
        T = D.'*spdiags(N,0,2,2)*D;
    else
        % a zero matrix - can set to scalar, will be taken as element-wise
        % addition (NOTE: this case should never be reached, added for safety)
      T = 0;
    end
end

% Find omega - accepts vectors!
function  omega = computeomega(Nleg,Mleg,eta)
    if ~isempty(eta)
        omega = 4./( 2*(Nleg+Mleg+eta)+1 )./( 2*(Nleg+Mleg+eta)+5 );
    else
        % Should never be reached, added for safety
        omega = 1;
    end
end

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
% END SCRIPT