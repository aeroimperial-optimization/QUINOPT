function INEQ = setBC(INEQ,BCEXPR)

%% SETBC.m Add boundary condition to integral inequality
%
% INEQ = SETBC(INEQ,BCEXPR) adds  the boundary condition BCEXPR=0 specified
% by expression EXPR to the integral inequality model INEQ. BCEXPR must be
% created using the dependent variables returned by QUIIMODEL. BCEXPR can be
% a vector of boundary conditions.
%
% See also BUILDINEQUALITYMODEL, DEPVAR

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    11/04/2015
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

if isempty(BCEXPR)
    % nothing to do
    return
end

% Get variables & base
BCEXPR = BCEXPR(:);
bcvars = getVariables(INEQ.BVAL,'stable');               % all boundary values variables
nbcvars = length(bcvars);

% exponents for rescaling if necessary
if INEQ.DOMAIN(1)~=-1 || INEQ.DOMAIN(2)~=1
    expon = cell2mat(arrayfun(@(x)(0:x),INEQ.MAXDER,'uniformoutput',0));    
    expon = transpose(vec(repmat(expon,2,1)));
end

for i=1:length(BCEXPR)
    
    [coeff,monom,ivars] = coefficients(BCEXPR(i));  % decompose bc
    [isbcvar,ind] = ismember(ivars,bcvars);         
    
    % Check the BC is ok
    if isZero(monom(1,:))
        % First monomial is always the constant 
        error('Boundary conditions are not homogeneous.')
        
    elseif ~all(isbcvar) || isa(coeff,'sdpvar') || isa(coeff,'legpoly')
        error(['BC must not depend on variables other than the boundary ',...
            'values of the dependent function.'])    
        
    elseif any( max(max(monom)) > 1 )
        error('Boundary conditions must be linear.')
        
    end
    
    % Find BC matrix & rescale to [-1,1] if necessary
    [I,J] = find(monom);
    coeff = coeff(I);
    if INEQ.DOMAIN(1)~=-1 || INEQ.DOMAIN(2)~=1
        d = ( 2/(INEQ.DOMAIN(2)-INEQ.DOMAIN(1))  ).^expon(ind(J));
        coeff = coeff.*d(:);
    end
    INEQ.BC = [INEQ.BC; sparse(1,ind(J),coeff,1,nbcvars)];
    
    
end

% Remove linearly dependent BCs
litol = 1e-10;
INEQ.BC = lirows(INEQ.BC,litol);


% ----------------------------------------------------------------------- %
%% NESTED FUNCTION - FIND DEPENDEND ROWS OF MATRIX

    function [Xsub,idx]=lirows(X,tol)
        % Extract a linearly independent set of rows of a given matrix X
        %
        %    [Xsub,idx]=lirows(X)
        %
        % inputs:
        %
        %  X: The given input matrix
        %  tol: A rank estimation tolerance. Default=1e-10
        %
        % outputs:
        %
        % Xsub: The extracted rows of X
        % idx:  The indices (into X) of the extracted rows
        
        % X has no non-zeros
        if ~nnz(X)
            Xsub=[]; idx=[];
            return
        end
        
        if nargin<2,
            tol=1e-10;
        end
        
        [Q,R,E] = qr(full(X'),0);   % full for compatibility with old MATLAB
        if ~isvector(R)
            diagr = abs(diag(R));
        else
            diagr = R(1);
        end
        
        r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
        deprows = size(X,1)-r;
        if deprows~=0
            fprintf('%i linearly dependent boundary conditions found and removed.\n',deprows);
        end
        idx=sort(E(1:r));
        Xsub=X(idx,:);
        
    end

% ----------------------------------------------------------------------- %
end

