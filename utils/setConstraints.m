function [CNSTR,SOS,LMIs,LINs,EQ,AuxVars,SoSParams] = setConstraints(x,Q,S,LMIs,slacks,AuxVars,EQ,opts)

%% SETCONSTRAINTS.m Set constraints for QuadIntIneq
%
% Set the constraints for the SDP relaxation of a quadratic integral
% inequality. Have SOS variables, linear matrix inequalities, linear
% inequalities, equalities and a list of auxiliary variables.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    05/05/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

CNSTR = [];

% ----------------------------------------------------------------------- %
% Infinite-dimensional part: SOS in [-1,1] or add to LMIs
% ----------------------------------------------------------------------- %
SOS = {};
SoSParams = [];
if degree(S,x)>0
    % Variables for S-procedure
    n = size(S,1);
    [N,C] = polySymmMatrix(x,opts.sosdeg,n);
    AuxVars{end+1,1} = C;
    if n>1;
        y = sdpvar(n,1);
    else
        y=1;
    end
    
    % Get SOS parameters
    SoSParams = findSosParameters(S,x);
    SoSParams = [SoSParams; C(:)];
    
    % Set SOS constraints
    CNSTR = [sos( y'*(S-N*(1-x^2))*y ); sos( y'*N*y )];
    SOS = {S-N*(1-x^2);N};
    
elseif isa(S,'sdpvar')
    % Add to LMIS
    LMIs{end+1,1} = S;
    
end


% ----------------------------------------------------------------------- %
% Linear matrix inequalities
% ----------------------------------------------------------------------- %
if isa(Q,'sdpvar');
    % Add Q to LMIs if variable, then set LMIs
    LMIs{end+1,1} = Q;
    for i=1:length(LMIs)
        CNSTR = [CNSTR; LMIs{i}>=0];
    end
    
elseif isnumeric(Q)
    % Test if it is positive definite
    E = eig((Q+Q')./2);
    if any(E<0)
        error('Infeasible problem (numeric LMI relaxation, not positive semidefinite).')
    end
    
    % Then set LMIs
    for i=1:length(LMIs)
        CNSTR = [CNSTR; LMIs{i}>=0];
    end
    
    % Then add numeric Q to list of LMIs
    LMIs{end+1,1} = Q;
end


% ----------------------------------------------------------------------- %
% Linear inequalities for slack variables
% ----------------------------------------------------------------------- %
LINs = [];
if ~isempty(slacks.t)
    LINs = [slacks.t - slacks.pcoef; slacks.t + slacks.pcoef];
    CNSTR = [CNSTR; LINs>=0];
end

% ----------------------------------------------------------------------- %
% Equality constraints from integration by parts
% ----------------------------------------------------------------------- %
if ~isempty(EQ)
    CNSTR = [CNSTR; EQ==0];
end

% END CODE
end