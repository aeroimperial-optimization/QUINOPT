function Z = any(X)

%% OVERLOADED: legpoly/any

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

[m,n] = size(X);
D = degree(X(:),[],'all');
ZeroDeg = find(D==0);
C = [X(ZeroDeg).coef];
if isa(C,'sdpvar')
    % Get sparsity pattern
    C = any(C);
end
Z = ones(m,n);
Z(ZeroDeg) = C~=0;
Z = reshape(Z,m,n);
Z = logical(Z);         % make logical