function INEQ = setIntegrand(INEQ,EXPR)

% SETINTEGRAND Set field F of structure INEQ
%
% INEQ = setIntegrand(INEQ,EXPR) rearranges the integrand of the integral 
%       inequality specified by EXPR as a matrix quadratic form and assigns
%       the computed matrices to the relevant fields INEQ.F.
%
% See also BUILDINEQUALITYMODEL

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    11/04/2015
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% Decompose expression in monomials
[coeff,monom,ivars] = coefficients(EXPR);

% Is EXPR homogeneous quadratic?
if any(sum(monom,2)~=2)
    error('EXPR must be a homogeneous quadratic polynomial of the dependent variables.')
end

% Extract used variables
allvars = getVariables([INEQ.DVAR; INEQ.BVAL],'stable');
ndvar = length(INEQ.DVAR);
Q = setQuadraticForm(coeff,monom,ivars,allvars);

% Rescale to [-1,1] if necessary
if INEQ.DOMAIN(1)~=-1 || INEQ.DOMAIN(2)~=1
    
    % Rescale independent variable inside Q (VERY SLOW FOR SDPVAR IF YALMIP HAS MANY INTERNALS)
    if isa(Q,'legpoly')
        Q = setDomain(Q,[-1,1]);
    elseif isa(Q,'sdpvar')
        z = ((INEQ.DOMAIN(2)-INEQ.DOMAIN(1))*INEQ.IVAR + INEQ.DOMAIN(2) + INEQ.DOMAIN(1))/2;
        Q = replace(Q,INEQ.IVAR,z);
    end
    
    % Rescale dependent variables & their derivatives
    d = cell2mat(arrayfun(@(x)(0:x),INEQ.MAXDER,'uniformoutput',0));    % list of exponents
    d = ( 2/(INEQ.DOMAIN(2)-INEQ.DOMAIN(1))  ).^d;
    d = [d, transpose(vec(repmat(d,2,1)))];
    Q = Q.*(d'*d);
    
end

% Find submatrices (Q is always upper triangular)
Fi = Q(1:ndvar,1:ndvar);
Fm = Q(1:ndvar,ndvar+1:end)';       % transpose since want in form BVAL*Fm*DVAR
Fb = Q(ndvar+1:end,ndvar+1:end);

% Set outputs if non-zero
if ~isZero(Fi); INEQ.F.Fi = Fi; end
if ~isZero(Fm); INEQ.F.Fm = Fm; end
if ~isZero(Fb);
    Fb = integrateBoundaryTerm(Fb,INEQ);        % must integrate beforehand!
    if isa(Fb,'legpoly'); Fb = sdpvar(Fb); end  % convert to sdpvar since no dependence on IVAR
    INEQ.F.Fb = Fb;
end

end

%%
% ----------------------------------------------------------------------- %
% NESTED FUNCTION
% ----------------------------------------------------------------------- %

function B = integrateBoundaryTerm(B,mod)
    % Integrate the boundary term to comply with internal model
    if isZero(B)
        return

    elseif degree(B,mod.IVAR)==0
        % B = 2.*sparse(double(B));     % BUG IF B IS SDPVAR?
        B = 2.*B;
        return

    elseif isa(B,'sdpvar')
        % Use YALMIP int function to integrate
        B = int(B,mod.IVAR,-1,1);
        return

    elseif isa(B,'legpoly')

        % Size and nonzero entries (B is square)
        n = size(B,2);
        [I,J] = find(any(B));

        % Loop to compute nonzero entries
        val = zeros(length(I),1);
        for k = 1:length(I)
            if degree(B(I(k),J(k))) == 0
                % Integrate a constant
                val(k) = 2*coefficients(B(I(k),J(k)));
            else
                % Integrate the coefficients
                c = coefficients(B(I(k),J(k)));

                % Integration matrix
                N = length(c);
                v1 = 1./(2*(0:N-1)'+1);
                v2 = [1;spalloc(N-1,1,0)];
                D = spdiags([v1 v2 -v1],-1:1,N+1,N);

                % Find integral: only need the coefficients in even position
                % (the others cancel out upon integration)
                c = D(2:2:end,:)*c(:);
                val(k) = 2*sum(c);

            end
        end

        % Set B
        B = sparse(I,J,val,n,n);
        return

    else
        error('Unknown class for INEQ.F.Fb.')

    end
end