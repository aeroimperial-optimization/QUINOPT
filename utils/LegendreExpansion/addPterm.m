function T = addPterm(T,Nleg,Mleg,pcoef,nnzIdx,dvar,ALPHA,BETA,DERORD)

%% addPterm.m
%
% [T,S,slk] = addPterm(T,S,slk,Nleg,Mleg,P,dvar,ALPHA,BETA,DERORD) computes
%   the relaxation of the R term of the Legendre expansion.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/10/2015
% Last Modified:    08/10/2015
% ----------------------------------------------------------------------- %

%% CODE

% ----------------------------------------------------------------------- %
% SETUP VARIABLES
% ----------------------------------------------------------------------- %
Ka = DERORD(dvar(1));   % maximum derivative of first variable
Kb = DERORD(dvar(2));   % maximum derivative of second variable

% Find the appropriate portion of the matrix T in which to add terms
TpartU = cumsum( Nleg+Mleg+1+2*DERORD ); % partition indices of T (upper limit)
TpartL = [1, TpartU(1:end-1)+1];         % partition indices of T (lower limit)
row = TpartL(dvar(1)):TpartU(dvar(1));   % row indices of block of T for the current variables
col = TpartL(dvar(2)):TpartU(dvar(2));   % col indices of block of T for the current variables

% ----------------------------------------------------------------------- %
% MATRIX OF INTEGRALS OF TRIPLE PRODUCTS OF LEGENDRE POLYNOMIALS
% ----------------------------------------------------------------------- %

if ALPHA==Ka && BETA==Kb
    % Both equal to the highest order derivative
    nMax = Nleg+Mleg+ALPHA; 
    mMax = Nleg+Mleg+BETA;
    LIMITS_ALPHA = [];
    LIMITS_BETA = [];
    
else
    % At most one highest order derivative (not both - already checked)
    nMax = Nleg+ALPHA;
    mMax = Nleg+BETA;
    
    % Set limits to construct integration matrices to handle case when
    % ALPHA or BETA are the highest order derivative
    if ALPHA==Ka
        LIMITS_ALPHA = [0, Nleg+Ka];
        LIMITS_BETA = [];
    elseif BETA==Kb
        LIMITS_ALPHA = [];
        LIMITS_BETA = [0, Nleg+Kb];
    else
        LIMITS_ALPHA = [];
        LIMITS_BETA = [];
    end
end

% Setup differentiation matrices
[Da,Ba] = legendreDiff(Nleg,Mleg,ALPHA,Ka,LIMITS_ALPHA);
[Db,Bb] = legendreDiff(Nleg,Mleg,BETA,Kb,LIMITS_BETA);

% Integrals of triple products - a cell array X
X = legendreTripleProduct(nnzIdx-1,0,nMax,0,mMax);

% ----------------------------------------------------------------------- %
% SET OUTPUT
% ----------------------------------------------------------------------- %
% Form the matrix representation of p(x)*d^ALPHA(u)*d^BETA(v)
Pmat = pcoef(1).*([Ba';Da']*X{1}*[Bb, Db]);
for j = 2:length(pcoef)
    % loop not entered if length(pcoef)<2!
    Pmat = Pmat + pcoef(j).*([Ba';Da']*X{j}*[Bb, Db]);
end

% Set output
T(row,col) = T(row,col) + Pmat;

end
%% END SCRIPT
