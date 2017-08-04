function [D,B] = legendreDiff(Nleg,Mleg,ALPHA,K,LIMITS,SYMM,options)

%% LEGENDREDIFF.m
%
% [D,B] = LEGENDREDIFF(Nleg,Mleg,ALPHA,K,LIMITS,SYMM,opts) returns the 
%       differentiation matrix that relates the coefficients of the derivative 
%       ALPHA of a function to those of the derivative K. The truncation of
%       the primitive function is assumed with coefficients from 0 to Nleg,
%       the derivative ALPHA is expanded with coefficients from 0 to Nleg+ALPHA
%       and the derivative K is expanded with coefficients from 0 to
%       Nleg+Mleg+K. When the optional input string "options" is set to 
%       'full', the function returns full matrices, otherwise the default 
%       outputs are sparse.
%
% NOTE: this version takes symmetry into account. Symmetry codes:
%       * 0: no symmetry
%       * 1: odd
%       * 2: even
%       Symmetry is for the zeroth derivative. The derivative ALPHA is symmetric
%       depending on the value of 
%       
%       ALPHA_DER_SYMM = (-1)^(ALPHA + SYMM);
%
%       If ALPHA_DER_SYMM = -1, then odd, if ALPHA_DER_SYMM = 1 then even. Set
%       ALPHA_DER_SYMM = 0 (no symmetry) if SYMM = 0.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    16/03/2015
% Last Modified:    24/04/2017
% ----------------------------------------------------------------------- %

%% CODE

% Set default if not specified
if nargin < 4
    error('Error in LegendreDiff.m: not enough input arguments.')
elseif nargin == 4
    LIMITS = [];
    SYMM = 0;
    options = 'sparse';
elseif nargin == 5
    SYMM = 0;
    options = 'sparse';
elseif nargin == 6
    options = 'sparse';
elseif nargin > 7
    error('Error in LegendreDiff.m: too many input arguments.')
end

% Get symmetry code
if SYMM==0
    ALPHA_DER_SYMM = 0;
    K_DER_SYMM = 0;
elseif SYMM==1 || SYMM==2
    ALPHA_DER_SYMM = (-1)^(ALPHA + SYMM);
    K_DER_SYMM = (-1)^(K + SYMM);
else
    error('Error in LegendreDiff2.m: invalid input SYMM (0, 1, or 2 allowed).')
end

% Construct matrices
if isempty(LIMITS)
    % ------------------------------------------------------------------- %
    % DIFFERENTIATION MATRICES WITH NO LIMITS
    % ------------------------------------------------------------------- %
    if ALPHA==K
        % No differentiation - just identity!
        D = speye(Nleg+Mleg+K+1);
        B = spalloc(Nleg+Mleg+K+1,K,0);
        
    else
        % Build differentiation recursively, use sparse matrices
        D = speye(Nleg+ALPHA+1);
        B = spalloc(Nleg+ALPHA+1,K,2^(K-ALPHA)-1);
        %S = sqrt( 2*(0:Nleg+ALPHA).'+1 )./sqrt(2);
        
        rt2 = sqrt(2);
        
        for i = 0:K-ALPHA-1
            
            eta = ALPHA+i;
            n = Nleg+eta+1;
            
            % matrix E_eta
            %             E = sparse(1,eta+1,1,n,K+2);
            %            B = B+D*E;
            
            % GF on 24 Apr 2017: the following is wrong - I fixed it
            % T = [sparse(Nleg+ALPHA+1,eta), ...
            %      D(:,1)./S, ...
            %      sparse(Nleg+ALPHA+1,K-1-eta)];
            T = [sparse(Nleg+ALPHA+1,eta), ...
                D(:,1).*rt2, ...
                sparse(Nleg+ALPHA+1,K-1-eta)];
            B = B + T;
            
            if n>=2
                % <-------
                % Old code - no scaling
                %C = [sparse([1 2],1,1,n,1), spdiags([v,-v],[-2 0],n,n)];
                % ------->
                
                % New code
                % Build rescaled matrix
                I = [1, 1:n, 2:n].';
                J = [1:n+1, 1:n-1].';
                v = [1, -1./(2*(1:n)+1), 1./(2*(0:n-2)+1)].';
                v = v.*( sqrt(2*(J-1)+1)./sqrt(2*(I-1)+1) );
                C = sparse(I,J,v,n,n+1);
                
            else
                % <-------
                % Old code - no scaling
                % C = [1, -1/3];
                % ------>
                
                % New code
                C = [1, -1/sqrt(3)];
                
            end
            
            D = D*C;
            
        end
        
        D = [D, spalloc(Nleg+ALPHA+1,Mleg,0)];
        
    end
    
    % Enforce symmetry if needed
    if ALPHA_DER_SYMM == -1 && K_DER_SYMM == -1
        % Set to zero even coefficients of both:
        % rows 1:2:end of B and D
        % cols 1:2:end of D
        D(1:2:end,:) = 0;
        D(:,1:2:end) = 0;
        B(1:2:end,:) = 0;
    elseif ALPHA_DER_SYMM == -1 && K_DER_SYMM == 1
        % Set to zero even coefficients of ALPHA and odd coeffs of K:
        % rows 1:2:end of B and D
        % cols 2:2:end of D
        D(1:2:end,:) = 0;
        D(:,2:2:end) = 0;
        B(1:2:end,:) = 0;
    elseif ALPHA_DER_SYMM == 1 && K_DER_SYMM == -1
        % Set to zero odd coefficients of ALPHA and even coeffs of K:
        % rows 2:2:end of B and D
        % cols 1:2:end of D
        D(2:2:end,:) = 0;
        D(:,1:2:end) = 0;
        B(2:2:end,:) = 0;
    elseif ALPHA_DER_SYMM == 1 && K_DER_SYMM == 1
        % Set to zero odd coefficients of ALPHA and odd coeffs of K:
        % rows 2:2:end of B and D
        % cols 2:2:end of D
        D(2:2:end,:) = 0;
        D(:,2:2:end) = 0;
        B(2:2:end,:) = 0;
    end
    
    
elseif (LIMITS(1)>=K-ALPHA)&&(LIMITS(2)<=Nleg+Mleg+ALPHA)
    % ------------------------------------------------------------------- %
    % DIFFERENTIATION MATRICES WITH LIMITS
    % ------------------------------------------------------------------- %
    % Also works for ALPHA=K - gives subset of identity matrix padded with zeros
    
    r = LIMITS(1);
    s = LIMITS(2);
    B = spalloc(s-r+1,K,0);     % no boundary terms!
    D = speye(s-r+1,s-r+1);
    
    for i = 0:K-ALPHA-1
        
        % matrix C^[r,s]
        
        % <-------
        % Old code - no scaling
        %         v = [1./(2*(r:s)'-1),  -1./(2*(r:s)'+3) ];
        %         C = spdiags(v,[0 2],s-r+1,s-r+3);
        % ------->
        
        % New code for rescaled matrix
        I = [1:s-r+1, 1:s-r+1].';
        J = [3:s-r+3, 1:s-r+1].';
        v = [-1./(2*(r:s)+3), 1./(2*(r:s)-1)].';
        v = v.*( sqrt(2*(J+r-2)+1)./sqrt(2*(I+r-1)+1) );
        C = sparse(I,J,v,s-r+1,s-r+3);
        
        D = D*C;
        r = r-1;    % update indices
        s = s+1;    % update indices
        
    end
    
    r = LIMITS(1);
    s = LIMITS(2);
    D = [spalloc(s-r+1,r+ALPHA-K,0), D, spalloc(s-r+1,Nleg+Mleg+ALPHA-s,0)];
    
    % Enforce symmetry if needed
    % Need to check if r is even or odd, and B already all zeros
    shift = 0;
    if rem(r,2)==1; shift=1; end
    
    if ALPHA_DER_SYMM == -1 && K_DER_SYMM == -1
        % Set to zero even coefficients of both:
        % rows 1+shift:2:end of D
        % cols 1:2:end of D
        D(1+shift:2:end,:) = 0;
        D(:,1:2:end) = 0;
    elseif ALPHA_DER_SYMM == -1 && K_DER_SYMM == 1
        % Set to zero even coefficients of ALPHA and odd coeffs of K:
        % rows 1+shift:2:end of D
        % cols 2:2:end of D
        D(1+shift:2:end,:) = 0;
        D(:,2:2:end) = 0;
    elseif ALPHA_DER_SYMM == 1 && K_DER_SYMM == -1
        % Set to zero odd coefficients of ALPHA and even coeffs of K:
        % rows 2-shift:2:end of D
        % cols 1:2:end of D
        D(2-shift:2:end,:) = 0;
        D(:,1:2:end) = 0;
    elseif ALPHA_DER_SYMM == 1 && K_DER_SYMM == 1
        % Set to zero odd coefficients of ALPHA and odd coeffs of K:
        % rows 2-shift:2:end of B and D
        % cols 2:2:end of D
        D(2-shift:2:end,:) = 0;
        D(:,2:2:end) = 0;
    end
    
    
else
    % ------------------------------------------------------------------- %
    % SOMETHING WRONG!
    % ------------------------------------------------------------------- %
    error(['Error in LegendreDiff.m: LIMITS should satisfy '...
        'LIMITS(1)>=K-ALPHA and LIMITS(2)<=Nleg+Mleg+ALPHA.'])
end

% Look at option to make D sparse or not
if strcmp(options,'full')
    D = full(D);
    B = full(B);
end



%% END CODE
end
