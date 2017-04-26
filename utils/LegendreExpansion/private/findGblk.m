function Gblk = findGblk(N,M,alpha,beta,k,rigorous,SYMM)

% FINDGBLK
%
% Part of QUINOPT
%
% Finds the matrix to expand the boundary values of the derivatives ALPHA to
% BETA of a function u expanded with N Legendre coefficients, in terms of the
% coefficients of the k-th derivative.
%
% NOTE: ALPHA,BETA < k!

if max(alpha,beta)>=k
    error('Inputs alpha and beta must be <= k-1')
elseif alpha>beta
    % swap around
    tmp = alpha;
    alpha = beta;
    beta = tmp;
end

if nargin<7
    SYMM=0;
end

% Get symmetry of k-th derivative
if SYMM==0
    K_DER_SYMM = 0;
elseif SYMM==1 || SYMM==2
    K_DER_SYMM = (-1)^(k + SYMM);
else
    error('Error in LegendreDiff2.m: invalid input SYMM (0, 1, or 2 allowed).')
end

% Initialise vectors for nonzero entries of G and their indices
dim = 2*(beta-alpha+1)+1+(2*k+1)*(beta-alpha)-2*sum(alpha+1:beta);
I = zeros(dim,1);
J = zeros(dim,1);
S = zeros(dim,1);

% Other variables: identity and zero vector
II = speye(k);
z = sparse(1,N+M+k+1);

% Build entries of G
count = 0;
rt2 = sqrt(2);

for der = alpha:min(beta,k-2)
    e = II(der+1,:);
    [D,B] = legendreDiff(N,M,der+1,k,[],SYMM);
    [Itmp,Jtmp,Stmp] = find([e, z; e+rt2.*B(1,:), rt2.*D(1,:)]);
    nInd = length(Itmp);                                        % how many nonzeros
    I(count+1:count+nInd) = Itmp + 2*(der-alpha);         % shift row indices
    J(count+1:count+nInd) = Jtmp;
    S(count+1:count+nInd) = Stmp;
    count = count+nInd;
end

if beta==k-1
    % Derivative k-1
    e = II(k,:);
    w = 1;
    if K_DER_SYMM == -1
        w = 0;
    end
    [Itmp,Jtmp,Stmp] = find([e, z; e, rt2.*w, z(2:end)]);
    nInd = length(Itmp);
    I(count+1:count+nInd) = Itmp + 2*(k-alpha) - 2;                 % shift row indices by 2*ALPHA-2
    J(count+1:count+nInd) = Jtmp;
    S(count+1:count+nInd) = Stmp;
end

% Remove potential zero indices
ind = ( (I.*J) > 0 );
I = I(ind);
J = J(ind);
S = S(ind);

% If not rigorous, set to zero entries with columns > N+1 (model fact that 
% the coefficient of monomials of order larger than N-k can be set to zero)
if ~rigorous
    ind = (J<=N+1);
    I = I(ind);
    J = J(ind);
    S = S(ind);
end

% Assemble G
Gblk = sparse(I,J,S,2*(beta-alpha+1),N+M+2*k+1);


end