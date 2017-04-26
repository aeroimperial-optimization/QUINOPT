function T = addPterm(T,Nleg,Mleg,pcoef,nnzIdx,dvar,ALPHA,BETA,DERORD,DVAR_SYMM,rigorous)

%% addPterm.m
%
% [T,S,slk] = addPterm(T,S,slk,Nleg,Mleg,P,dvar,ALPHA,BETA,DERORD,rigorous) 
%   computes the relaxation of the P term of the Legendre expansion. If
%
%       rigorous = 0
%
%   then also set entries of the bottom-right corner of T to zero, to model the
%   fact that the expansion of u with N coefficients only is finite dimensional 
%   and the Legendre coefficients of the k-th order derivative from N-k+1 to
%   N+Mleg+k should be set to zero. That is, keep only the following N+1
%   variables for each dependent variable:
%
%   - k boundary values
%   - N-k+1 Legendre coefficients (with indices from 0 to N-k)


% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/10/2015
% Last Modified:    03/06/2016
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
[Da,Ba] = legendreDiff(Nleg,Mleg,ALPHA,Ka,LIMITS_ALPHA,DVAR_SYMM(dvar(1)));
[Db,Bb] = legendreDiff(Nleg,Mleg,BETA,Kb,LIMITS_BETA,DVAR_SYMM(dvar(2)));

% Integrals of triple products - a cell array X
X = legendreTripleProduct(nnzIdx-1,0,nMax,0,mMax);

% ----------------------------------------------------------------------- %
% SET OUTPUT
% ----------------------------------------------------------------------- %
% Form the matrix representation of p(x)*d^ALPHA(u)*d^BETA(v)
Ma = [Ba.'; Da.'];
Mb = [Bb, Db];
Pmat = pcoef(1).*(Ma*X{1}*Mb);
for j = 2:length(pcoef)
    % loop not entered if length(pcoef)<2!
    Pmat = Pmat + pcoef(j).*(Ma*X{j}*Mb);
end

% Set entries to zero if not rigorous
% ADDED 
if ~rigorous
    Pmat(Nleg+2:end,:) = 0;     % set last rows to zero
    Pmat(:,Nleg+2:end) = 0;     % set last columns to zero
end

% Assign output
if isa(T,'sdpvar') || isa(Pmat,'sdpvar')
    T = sdpvarAddInPlace(T,Pmat,row,col);
else
    T(row,col) = T(row,col) + Pmat;
end

end
%% END SCRIPT
