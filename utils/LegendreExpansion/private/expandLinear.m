function L = expandLinear(INEQ,Nleg,Mleg,P,opts)

% EXPANDLINEAR
%
% Expand the linear part of the integrand in the integral inequality. Returns a
% row vector for convenience.

% Extract useful variables
dimint = INEQ.dimint;
dimbnd = INEQ.dimbnd;

% Do we need to expand at all?
if ( isempty(INEQ.L.Lb) || isZero(INEQ.L.Lb) ) && ( isempty(INEQ.L.Li) || isZero(INEQ.L.Li))
    L = sparse(1,dimint+dimbnd);
    return
    
else
    
    % Expand boundary term - distinction between rigorous
    % and non-rigorous expansion already taken into account by matrix P.
    if ~isempty(INEQ.L.Lb) && ~isZero(INEQ.L.Lb);
        L = INEQ.L.Lb.'*P;
    end
    
    
    % Expand mixed term - distinction between rigorous and non-rigorous expansion
    % taken into account by matrix P and when building Qmix by setting entries to 0.
    % However, still need to fix the outer approximation if not rigorous: expand
    % with more coefficient than in "expandLi.m".
    if ~isempty(INEQ.L.Li) && ~isZero(INEQ.L.Li);
        Lint = expandLi(Nleg,Mleg,INEQ.L.Li,INEQ.IVAR,INEQ.DERORD,opts.rigorous);
        if opts.rigorous
            Lint = [Lint, sparse(1,dimbnd)];
        else
            Lint = fixOuterApproxLi(Lint,INEQ,Nleg,Mleg);
        end
        L = L + Lint;
    end
    
end