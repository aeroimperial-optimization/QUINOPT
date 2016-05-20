function varargout = quinopt(EXPR,BC,OBJ,CNSTR,PARAMETERS,N,OPTIONS)

%% QUINOPT.m Solve optimization problem with quadratic integral inequalities
%
% SOL = QUINOPT(EXPR) tests whether the q homogeneous quadratic integral 
%     inequalities
%
%     \int_{a}^{b} EXPR(i) dx >=0, i = 1, ..., q
%
%     are feasible using a finite dimensional relaxation based on semidefinite
%     programming. As shown in the expression above, EXPR is a vector whose
%     i-th entry specifies the integrand of the i-th integral inequality.
%     Each entry of EXPR must be a quadratic homogeneous polynomial of the
%     variables returned by the commands <a href="matlab:help('indvar')">indvar</a> and <a href="matlab:help('depvar')">depvar</a> variables. 
%     Examples can be found in the folder "QUINOA_vX.Y/resources", where X.Y
%     is your version number. The output SOL is a structure with the
%     following fields:
%
%     - setupTime: the time taken to set up the problem
%     - solutionTime: the time taken by YALMIP to solve the problem
%     - problem: code of problem encountered during setup. Values are:
%                * 0: no problem
%                * 1: ill-posed inequality
%                * 2: infeasible relaxation 
%     - FeasCode: code for the feasibility of the solution returned by YALMIP. 
%                 See <a href="matlab:help('quinoptFeasCode')">quinoptFeasCode</a> for a complete list.
%     - YALMIP: the solution structure returned by YALMIP. See <a href="matlab:help('optimize')">optimize</a>
%               and <a href="matlab:help('yalmip/modules/sos/solvesos')">solvesos</a> for more details.
%
% SOL = QUINOPT(EXPR,BC) determines whether the integral inequalities
%     specified by EXPR are feasible over the set defined by the homogeneous
%     boundary conditions specified by the vector BC. Specifically, BC is
%     interpreted as the list of boundary conditions BC(1)=0, ..., BC(end)=0.
%     Like EXPR, BC must be created using the variables returned  by the 
%     commands <a href="matlab:help('indvar')">indvar</a> and <a href="matlab:help('depvar')">depvar</a>.
%
% SOL = QUINOPT(EXPR,BC,OBJ) optimizes the objective function OBJ
%     constrained by the integral inequalities specified by EXPR.
%
% SOL = QUINOPT(EXPR,BC,OBJ,CNSTR) and SOL = QUINOPT(EXPR,BC,OBJ,CNSTR,PARAMETERS)
%     optimize the objective function OBJ subjet to the integral
%     inequalities specified by EXPR and BC, and the additional constraints
%     given by CNSTR. If CNSTR contains sum-of-square constraints, then
%     the parameters in the polynomial expressions MUST be specified in
%     the input vector PARAMETERS. See <a href="matlab:help('@sdpvar/sos')">sos</a>  and  <a href="matlab:help('yalmip/modules/sos/solvesos')">solvesos</a> for more details 
%     on specifying sum-of-squares constraints with YALMIP syntax.
%
% SOL = QUINOPT(EXPR,BC,OBJ,CNSTR,PARAMETERS,N) also specifies the number of
%     Legendre coefficients to use in the relaxation of the integral
%     inequality.
%
% SOL = QUINOPT(EXPR,BC,OBJ,CNSTR,PARAMETERS,N,OPTIONS) overrides the
%     default options. OPTIONS is a structure containing the following
%     fields:
%
%     - YALMIP: a substructure containing the options for YALMIP, set
%               with YALMIP's command <a href="matlab:help('sdpsettings')">sdpsettings</a>.
%     - rigorous: if set to 'true', the semidefinite relaxation estimates
%               the contribution of the infinitly many Legendre modes of the 
%               dependent variables of degree larger than N. If set to 
%               'false', the Legendre series expansions are simply truncated.
%     - BCprojectorBasis: string specifying which basis to use for the 
%               projection on the boundary conditions. If set to 'rref', 
%               use a "rational" basis. If set to 'orth', use an 
%               orthonormal basis. The orthonormal basis may be preferable
%               numerically, but it may destroy sparsity of the data.
%     - sosdeg: the degree of the sum-of-squares polynomials used in the
%               S-procedure to localize SOS constraints from the integral
%               inequality to the integration domain.
%
% [SOL,CNSTR,DATA] = QUINOPT(...) also returns the YALMIP constraint object
%      CNSTR used to solve the optimization problem, and a structure DATA
%      containing all raw variables used to set up the constraints in
%      CNSTR.
%
% Unused inputs can be left empty; for example, QUINOPT(EXPR,BC,[],CNSTR)
% determines whether the integral inequality and the constraints in CNSTR
% are feasible without optimizing an objective function.
%
% See also INDVAR, DEPVAR, SETQUADINTINEQ, OPTIMIZE, @SDPVAR/VALUE, SDPSETTINGS, 
%          CLEARMODEL, QUINOPTFEASCODE

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    05/05/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% Set empty arguments
if nargin > 7; error('Too many inputs.'); end
if nargin < 7; OPTIONS = struct([]); end
if nargin < 6; N = []; end
if nargin < 5; PARAMETERS = []; end
if nargin < 4; CNSTR = []; end
if nargin < 3; OBJ = []; end
if nargin < 2; error('Not enough inputs.'); end

% Check on inputs
if ~isempty(EXPR) && ~isa(EXPR,'dvarpoly')
    error('Class of EXPR must be "dvarpoly". Please use the commands INDVAR and DEPVAR to create your problem variables.')
elseif ~isempty(BC) && ~isa(BC,'dvarpoly')
    error('Class of BC must be "dvarpoly". Please use the commands INDVAR and DEPVAR to create your problem variables.')
elseif ~isempty(OBJ) && numel(OBJ)>1
    error('Objective function should be a scalar expression.')
end
    
%Set options
OPTIONS = setQuadIntIneqOptions(OPTIONS);

% Setup integral inequality constraints
time = tic;
for i = length(EXPR):-1:1

    [QIICNSTR,DATA(i),FLAG] = setQuadIntIneq(EXPR(i),BC,N,OPTIONS);
    SOL.setupTime = toc(time);
    SOL.solutionTime = [];
    SOL.problem = FLAG;

    if FLAG==0
        % No problem, add constraints to list
        CNSTR = [QIICNSTR; CNSTR];
        PARAMETERS = [PARAMETERS; DATA(i).SumOfSquaresParameters];
    else
        % Problem infeasible (see warning message displayed for reason)
        SOL.FeasCode = 1;
        SOL.YALMIP = [];
        return
    end

end

% Find if have a SOS problem
try 
    issos = any(is(CNSTR,'sos'));
catch
    issos = 0;
end

% Solve if no problem during setup
if FLAG==0 && issos
    time = tic;
    [yalmipsol,m,Q,res,everything] = solvesos(CNSTR,OBJ,OPTIONS.YALMIP,PARAMETERS);
    SOL.solutionTime = toc(time);
    SOL.FeasCode = yalmipsol.problem;
    SOL.YALMIP = yalmipsol;
    SOL.YALMIP.monomials = m;
    SOL.YALMIP.sosDecompositionMatrices = Q;
    SOL.YALMIP.residuals = res;
    SOL.YALMIP.everything = everything;
    
else
    time = tic;
    yalmipsol = optimize(CNSTR,OBJ,OPTIONS.YALMIP);
    SOL.solutionTime = toc(time);
    SOL.FeasCode = yalmipsol.problem;
    SOL.YALMIP = yalmipsol;
    
end

% Set outputs
if nargout > 0; varargout{1} = SOL; end;
if nargout > 1; varargout{2} = CNSTR; end;
if nargout > 2; varargout{3} = DATA; end;
if nargout > 3; error('You have requested too many outputs.'); end;

% END CODE
end