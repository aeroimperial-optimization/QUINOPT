function [CNSTR,DATA,FLAG] = setConstraints(x,Q,S,LMIs,slacks,AuxVars,EQ,FLAG,opts)

%% SETCONSTRAINTS.m Set constraints for QuadIntIneq
%
% Set the constraints for the SDP relaxation of a quadratic integral
% inequality. Outputs:
%
% CNSTR: A structure containing the YALMIP constraints to pass to optimize. An
%        empty structure is returned if one of the constraint is detected to be
%        infeasible. A list of infeasible constraints is given in
%        DATA.InfeasibleConstraints (see below)
%
%  DATA: A structure containing the raw data used to set up the CNSTR structure.
%        The fields are:
%       - SumOfSquares          : a list of polynomials to be sum-or-squares
%       - SumOfSquaresParameters: the optimization variables in the polynomials
%                                 listed in DATA.SumOfSquares
%       - MatrixInequalities    : a list of matrices to be positive semidefinite
%       - LinearInequalities    : a list of expressions to be nonnegative
%       - Equalities            : a list of expressions to be equal to zero
%       - slacks                : a list of the slack variables used to handle
%                                 absolute values
%       - auxiliaryVariables    : a list of auxiliary variables introduced by
%                                 QUINOPT in addition to the original problem
%                                 variables
%       - InfeasibleConstraints : a list of variables that should be positive
%                                 semidefinite but are not, proving that the 
%                                 problem is trivially infeasible.
%
%  FLAG: An output flag: 1 if all good, 2 if at least one infeasible constraint 
%        was detected.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    05/05/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %


% ----------------------------------------------------------------------- %
% Useful stuff
% ----------------------------------------------------------------------- %
CNSTR = [];                 % the constraint structure
INFCNSTR = {};              % list of infeasible constraints
warnstr = sprintf(['\nWARNING: The current relaxation is infeasible. Try increasing the\n',...
             'number of Legendre coefficients used in the expansion using the  \n',...
             'input N. If the problem persists, your problem might be infeasible\n',...
             'although a certificate of infeasibility is not available.\n\n']);
         
         
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
    % Add to LMIs
    LMIs{end+1,1} = S;
    
elseif isnumeric(S)
    % Check for feasibility, display warning if not feasible and set 'solve'
    % option to 'false'
    if any(eig(S)<opts.psdtol)
        INFCNSTR{end+1,1} = S;
        fprintf(warnstr);
        FLAG = 2;
    end
    
    % Add to LMIs
    LMIs{end+1,1} = S;
    
end


% ----------------------------------------------------------------------- %
% Finite dimensional part (the large LMI)
% ----------------------------------------------------------------------- %
if isa(Q,'sdpvar');
    % Add Q to LMIs if variable, then set LMIs
    LMIs{end+1,1} = Q;
    
elseif isnumeric(Q)
    % Test if it is positive definite
    E = eig((Q+Q')./2);
    if any(E<opts.psdtol)
        INFCNSTR{end+1,1} = Q;
        fprintf(warnstr);
        FLAG = 2;
    end
        
    % Then add numeric Q to list of LMIs
    LMIs{end+1,1} = Q;
end


% ----------------------------------------------------------------------- %
% Linear matrix inequalities
% ----------------------------------------------------------------------- %
% Do not set the numeric LMIs - already tested!
% Only set constraints if no numerical LMIs are infeasible
if FLAG==0
    for i=1:length(LMIs)
        if isa(LMIs{i},'sdpvar')
            CNSTR = [CNSTR; LMIs{i}>=0];
        end
    end
end

% ----------------------------------------------------------------------- %
% Linear inequalities for slack variables
% ----------------------------------------------------------------------- %
% Only set constraints if FLAG==0 (no problems detected, worth trying to solve)
LINs = [];
if ~isempty(slacks.t)
    LINs = [slacks.t - slacks.pcoef; slacks.t + slacks.pcoef];
    if FLAG==0
        CNSTR = [CNSTR; LINs>=0];
    end
end


% ----------------------------------------------------------------------- %
% Equality constraints from integration by parts
% ----------------------------------------------------------------------- %
% Only set constraints if FLAG==0 (no problems detected, worth trying to solve)
if ~isempty(EQ) && FLAG==0
    CNSTR = [CNSTR; EQ==0];
end

% ----------------------------------------------------------------------- %
% Set DATA structure
% ----------------------------------------------------------------------- %
DATA.SumOfSquares = SOS;
DATA.SumOfSquaresParameters = SoSParams;
DATA.MatrixInequalities = LMIs;
DATA.LinearInequalities = LINs;
DATA.Equalities = EQ;
DATA.slacks = slacks;
DATA.auxiliaryVariables = AuxVars;
DATA.InfeasibleConstraints = INFCNSTR;

% END CODE
end