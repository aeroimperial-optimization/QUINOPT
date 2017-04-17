function INEQ = setInequalityModel(EXPR,BC)

%% SETINEQUALITYMODEL.m Build inequality model
%
% INEQ = SETINEQUALITYMODEL(EXPR,BC) constructs the internal model used
%       to represent the integral inequality specified by the integrand
%       EXPR and the boudnary conditions BC.
%
% See also ADDBC, SETINTEGRAND

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    11/04/2015
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% Check inputs
if isempty(EXPR)
    INEQ = struct([]);
    return
elseif ~isa(EXPR,'dvarpoly') || ~isscalar(EXPR)
    error('EXPR must be a 1-by-1 polynomial expression of the dependent variables (class dvarpoly).')
elseif ~isempty(BC) && (~isa(BC,'dvarpoly') || ~isvector(BC))
    error('BC must be a vector of the boundary values of the dependent variables (class dvarpoly).')
end

% Get ID of dependent variables in expression and in internal model
ivars = getVariables([EXPR;BC(:)],'stable');
mod = qiimodel('query');
if isempty(mod.DEPVARMODEL)
    error('Internal model is empty. Please setup your problem using the commands "indvar" and "depvar".')
else
    dvarID = getVariables(mod.DEPVARMODEL.DVAR,'stable');
    bvalID = getVariables(mod.DEPVARMODEL.BVAL,'stable');
end

% Find if dependent variables in internal model is used and internal ID of
% used variables.
varsets = mat2cell([dvarID,bvalID],1,[mod.DEPVARMODEL.MAXDER+1, 2*mod.DEPVARMODEL.MAXDER+2]);
varsets = reshape(varsets,length(mod.DEPVARMODEL.MAXDER),2);
isUsed = any(cellfun(@(x)any(ismember(ivars,x)),varsets),2);
dvarNewID = [varsets{isUsed,1}];        % identifiers of used dependent variables
bvalNewID = [varsets{isUsed,2}];        % identifiers of used boundary variables

% Check that only one dependent variable was used
IVAR = mod.DEPVARMODEL.IVAR(isUsed);
if length(depends(IVAR)) > 1
    error(['Your problem has %i independent variables, but QuadIntIneq ',...
           'only supports problems with one independent variable.'],length(depends(IVAR)));
end

% Set inequality model
INEQ.IVAR = IVAR(1);
INEQ.DOMAIN = mod.INDVARMODEL.DOMAIN(mod.INDVARMODEL.IVARID==depends(IVAR),:);
INEQ.MAXDER = mod.DEPVARMODEL.MAXDER(isUsed);
INEQ.DERORD = INEQ.MAXDER;
INEQ.DVAR = mod.DEPVARMODEL.DVAR(ismember(dvarID,dvarNewID));
INEQ.BVAL = mod.DEPVARMODEL.BVAL(ismember(bvalID,bvalNewID));
INEQ.F.Fi = [];
INEQ.F.Fm = [];
INEQ.F.Fb = [];
INEQ.L.Li = [];
INEQ.L.Lb = [];
INEQ.C = [];
INEQ.BC = [];

% Set fields F and BC using dedicated functions
INEQ = setIntegrand(INEQ,EXPR);            
INEQ = setBC(INEQ,BC);  

%% END CODE
end