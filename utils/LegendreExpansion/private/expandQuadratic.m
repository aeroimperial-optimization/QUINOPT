function [Q,S,slacks,MatrixInequalities,AuxVars] = ...
    expandQuadratic(INEQ,Nleg,Mleg,isFi,isFm,isFb,P,opts)

% EXPANDQUADRATIC
%
% Expand the quadtic part of the integrand in the integral inequality.

% Extract useful variables
dimint = INEQ.dimint;
dimbnd = INEQ.dimbnd;

% Expand integral term - distinction between rigorous and non-rigorous expansion 
% taken into account when building Q by setting entries to 0.
Q(dimint,dimint) = INEQ.IVAR;  % initialize, fake dependence on IVAR
if isFi
    [Q,S,slacks,MatrixInequalities,AuxVars] = ...
        expandFi(Q,Nleg,Mleg,INEQ.F.Fi,INEQ.IVAR,INEQ.DERORD,opts);
end
Q = replace(Q,INEQ.IVAR,0); % remove fake dependence on IVAR
if opts.rigorous
    Q = [Q, sparse(dimint,dimbnd); sparse(dimbnd,dimint), sparse(dimbnd,dimbnd)];
else
    [Q,Fix] = fixOuterApproxQ(Q,INEQ,Nleg,Mleg);
end


% Expand boundary term (need to make symmetric) - distinction between rigorous
% and non-rigorous expansion already taken into account by matrix P.
if isFb
    Qbnd = P'*( (INEQ.F.Fb + INEQ.F.Fb.')./2 )*P;
    Q = Q + Qbnd;
end


% Expand mixed term - distinction between rigorous and non-rigorous expansion 
% taken into account by matrix P and when building Qmix by setting entries to 0.
if isFm
    Qmix = expandFm(Nleg,Mleg,INEQ.F.Fm,INEQ.IVAR,INEQ.DERORD,opts.rigorous);
    Qmix = P'*Qmix;
    if opts.rigorous
        Qmix = [Qmix, sparse(dimint+dimbnd,dimbnd)];
    else
        Qmix = Qmix*Fix;
    end
    Q = Q + 0.5.*(Qmix+Qmix.');
end