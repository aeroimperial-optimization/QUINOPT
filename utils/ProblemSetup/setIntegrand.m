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
mon_deg = sum(monom,2);
if any(mon_deg>2)
    error('EXPR must be a quadratic polynomial of the dependent variables.')
end

% Extract used variables
allvars = getVariables([INEQ.DVAR; INEQ.BVAL],'stable');
ntotvars = length(allvars);
ndvar = length(INEQ.DVAR);

% Set homogeneous quadratic term
quadratic = (mon_deg==2);
if any(quadratic)
    Q = setQuadraticForm(coeff(quadratic),monom(quadratic,:),ivars,allvars);
else
    Q = zeros(ntotvars);
end

% Set linear term
linear = (mon_deg==1);
if any(linear)
    L = setLinearForm(coeff(linear),monom(linear,:),ivars,allvars);
else
    L = zeros(ntotvars,1);
end

% Set constant
const = (mon_deg==0);
if any(const)
    C = coeff(const);
else
    C = 0;
end

% Rescale to [-1,1] if necessary
if INEQ.DOMAIN(1)~=-1 || INEQ.DOMAIN(2)~=1
    
    % Rescale independent variable inside Q (VERY SLOW FOR SDPVAR IF YALMIP HAS MANY INTERNALS)
    if isa(Q,'legpoly')
        Q = setDomain(Q,[-1,1]);
    elseif isa(Q,'sdpvar')
        z = ((INEQ.DOMAIN(2)-INEQ.DOMAIN(1))*INEQ.IVAR + INEQ.DOMAIN(2) + INEQ.DOMAIN(1))/2;
        Q = replace(Q,INEQ.IVAR,z);
    end
    
    % Rescale independent variable inside L (VERY SLOW FOR SDPVAR IF YALMIP HAS MANY INTERNALS)
    if isa(L,'legpoly')
        L = setDomain(L,[-1,1]);
    elseif isa(L,'sdpvar')
        z = ((INEQ.DOMAIN(2)-INEQ.DOMAIN(1))*INEQ.IVAR + INEQ.DOMAIN(2) + INEQ.DOMAIN(1))/2;
        L = replace(L,INEQ.IVAR,z);
    end
    
    % Rescale independent variable inside C (VERY SLOW FOR SDPVAR IF YALMIP HAS MANY INTERNALS)
    if isa(C,'legpoly')
        C = setDomain(C,[-1,1]);
    elseif isa(C,'sdpvar')
        z = ((INEQ.DOMAIN(2)-INEQ.DOMAIN(1))*INEQ.IVAR + INEQ.DOMAIN(2) + INEQ.DOMAIN(1))/2;
        C = replace(C,INEQ.IVAR,z);
    end
    
    % Rescale dependent variables & their derivatives
    d = cell2mat(arrayfun(@(x)(0:x),INEQ.MAXDER,'uniformoutput',0));    % list of exponents
    d = ( 2/(INEQ.DOMAIN(2)-INEQ.DOMAIN(1))  ).^d;
    d = [d, transpose(vec(repmat(d,2,1)))];
    Q = Q.*(d.'*d);
    L = L.*d(:);
    
end

% Find submatrices of Q (Q is always upper triangular)
Fi = Q(1:ndvar,1:ndvar);
Fm = Q(1:ndvar,ndvar+1:end)';       % transpose since want in form BVAL*Fm*DVAR
Fb = Q(ndvar+1:end,ndvar+1:end);

% Find subvectors of L 
Li = L(1:ndvar);
Lb = L(ndvar+1:end);

% Set outputs if non-zero
if ~isZero(Fi); INEQ.F.Fi = Fi; end
if ~isZero(Fm); INEQ.F.Fm = Fm; end
if ~isZero(Fb)
    Fb = integrateBoundaryTerm(Fb,INEQ);        % must integrate beforehand!
    if isa(Fb,'legpoly'); Fb = sdpvar(Fb); end  % convert to sdpvar since no dependence on IVAR
    INEQ.F.Fb = Fb;
end

if ~isZero(Li); INEQ.L.Li = Li; end
if ~isZero(Lb)
    Lb = integrateBoundaryTerm(Lb,INEQ);        % must integrate beforehand!
    if isa(Lb,'legpoly'); Lb = sdpvar(Lb); end  % convert to sdpvar since no dependence on IVAR
    INEQ.L.Lb = Lb;
end

if ~isZero(C)
    C = integrateBoundaryTerm(C,INEQ);        % must integrate beforehand!
    if isa(C,'legpoly'); C = sdpvar(C); end  % convert to sdpvar since no dependence on IVAR
    INEQ.C = C;
elseif isZero(C)
    INEQ.C = 0;
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
        [m,n] = size(B);
        [I,J] = find(any(B));
    
        % Loop to compute nonzero entries
        % val = zeros(length(I),1);
        for k = length(I):-1:1
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
        B = sparse(I,J,val,m,n);
        return

    else
        error('setIntegrand.integrateBoundaryTerm(B,mod): unknown class for B.')

    end
end