function [Q,H] = fixOuterApproxLi(Q,INEQ,N,M)

% FIXOUTERAPPROXL
%
% Fix L for outer aproximation by expanding everything in terms of derivatives
% of order specified by INEQ.MAXDER (plus one to enhance sparsity)


% Loop over dependent variables
H = [];
for i = 1:INEQ.ndvars
    Hblk = findHblk(INEQ.DERORD(i),INEQ.MAXDER(i),N,M);
    H = spblkdiag(H,Hblk);
end

% Fix Q
Q = Q*H;

end


%% Nested function

function Hblk = findHblk(k,l,N,M)

% Find matrix to expand variables of k-th derivative to l-th derivative
[D,B] = legendreDiff(N,M,k,l+1); 
E = spdiags(ones(k,1),0,k,l+1);
Hblk = [E, sparse(k,size(D,2)); B, D];

% Since not rigorous, set to zero entries with columns > N+1 (model fact that 
% the coefficient of monomials of order larger than N-l-1 can be set to zero)
Hblk(:,N+2:end)=0;

end