function [D,B] = legendreDiff(Nleg,Mleg,ALPHA,K,LIMITS,options)

%% LEGENDREDIFF.m
%
% [D,B] = LEGENDREDIFF(Nleg,Mleg,ALPHA,K,LIMITS,opts) returns the 
%       differentiation matrix that relates the coefficients of the derivative 
%       ALPHA of a function to those of the derivative K. The truncation of
%       the primitive function is assumed with coefficients from 0 to Nleg,
%       the derivative ALPHA is expanded with coefficients from 0 to Nleg+ALPHA
%       and the derivative K is expanded with coefficients from 0 to
%       Nleg+Mleg+K. When the optional input string "options" is set to 
%       'full', the function returns full matrices, otherwise the default 
%       outputs are sparse.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    16/03/2015
% Last Modified:    05/04/2016
% ----------------------------------------------------------------------- %

%% CODE

% Set default if not specified
if nargin < 4
    error('Error in LegendreDiff.m: not enough input arguments.')
elseif nargin == 4
    LIMITS = [];
    options = 'sparse';
elseif nargin == 5
    options = 'sparse';
elseif nargin > 6
    error('Error in LegendreDiff.m: too many input arguments.')
end

% Construct matrices
if isempty(LIMITS);
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
        
        for i = 0:K-ALPHA-1
            
            eta = ALPHA+i;
            n = Nleg+eta+1;
            
            % matrix E_eta
            %             E = sparse(1,eta+1,1,n,K+2);
            %             B = B+D*E;
            B = B + [spalloc(Nleg+ALPHA+1,eta,0),D(:,1),spalloc(Nleg+ALPHA+1,K-1-eta,0)];
            
            % matrix C_eta
            v = 1./(2*(0:Nleg+eta)'+3);
            if n>=2
                C = [sparse([1 2],1,1,n,1), spdiags([v,-v],[-2 0],n,n)];
            else
                C = [1, -1/3];
            end
            
            D = D*C;
            
        end
        
        D = [D, spalloc(Nleg+ALPHA+1,Mleg,0)];
        
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
        v = [1./(2*(r:s)'-1),  -1./(2*(r:s)'+3) ];
        C = spdiags(v,[0 2],s-r+1,s-r+3);
        
        D = D*C;
        r = r-1;    % update indices
        s = s+1;    % update indices
        
    end
    
    r = LIMITS(1);
    s = LIMITS(2);
    D = [spalloc(s-r+1,r+ALPHA-K,0), D, spalloc(s-r+1,Nleg+Mleg+ALPHA-s,0)];
    
    
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
