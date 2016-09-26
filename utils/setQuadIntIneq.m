function varargout = setQuadIntIneq(EXPR,BC,N,OPTIONS)

%% SETQUADINTINEQ.m Create constraint for quadratic integral inequality
%
% CNSTR = SETQUADINTINEQ(EXPR,BC) returns a YALMIP constraint object
%       representing a semidefinite relaxation of a homogeneous quadratic
%       integral inequality constraint. The integrand of the inequality is 
%       specified by EXPR, and the boundary conditions on the dependent 
%       variables are specified by BC. Specifically, the inputs are
%
%       - EXPR: a symbolic expression defining the integrand of the
%               inequality constraint. EXPR must be created using the variables
%               returned by the functions <a href="matlab:help('indvar')">indvar</a> and <a href="matlab:help('depvar')">depvar</a>.
%
%       - BC: a vector of symbolic expression that define the boundary
%             conditions. These are interpreted as BC(1)=0, ..., BC(end)=0.
%             BC must be created using the variables returned by <a href="matlab:help('depvar')">depvar</a>.
%       
%       - N (optional): The number of Legendre coefficients used to formulate
%                       the SDP relaxation of the integral inequality.
%                       Larger values of N result in better relaxations,
%                       but also give larger optimization problems.
%
%       - OPTIONS: a structure containing options for SETQUADINTINEQ. Allowed
%                  fields are:
%
%               - YALMIP: a substructure containing the options for YALMIP,
%                   set with YALMIP's command <a href="matlab:help('sdpsettings')">sdpsettings</a>.
%               - rigorous: if set to 'true', the semidefinite relaxation 
%                   estimates the contribution of the infinitly many Legendre
%                   modes of the dependent variables of degree larger than 
%                   N. If set to 'false', the Legendre series expansions are
%                   simply truncated.
%               - BCprojectorBasis: string specifying which basis to use for
%                   the projection on the boundary conditions. If set to 
%                   'rref', use a "rational" basis. If set to 'orth', use an 
%                   orthonormal basis. The orthonormal basis may be preferable
%                   numerically, but it may destroy sparsity of the data.
%               - sosdeg: the degree of the sum-of-squares polynomials used
%                   in the S-procedure to localize SOS constraints arising
%                   from the relaxation of the integral inequality to the 
%                   domain of integration.
%               - psdtol: (default -1e-8) a tolerance to establish if a negative
%                   eigenvalue of a matrix should be considered negative, or if
%                   can be considered as a numerical zero with roundoff error.
%
% [CNSTR,DATA,FLAG] = SETQUADINTINEQ(EXPR,BC) also returns the raw data used to
%       set up the constraint object CNSTR and a FLAG to determine if and
%       which problem occurred. Values for FLAG are:
%       - 0: no problem
%       - 1: ill-posed inequality
%       - 2: infeasible relaxation 
%
% See also  INDVAR, DEPVAR, QUINOPT, CLEARMODEL

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    05/05/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

%% CODE

% Initialize output
varargout = cell(nargout,1);

% ----------------------------------------------------------------------- %
% Set default inputs
% ----------------------------------------------------------------------- %
if nargin == 0
        error('Not enough input arguments.')
elseif nargin >= 1 && isempty(EXPR)
       if nargout > 2; varargout{3} = 0; end       % FLAG
       return
elseif nargin == 1
        BC = [];
        N = [];
        OPTIONS = [];
elseif nargin == 2
        N = [];
        OPTIONS = [];
elseif nargin == 3
        OPTIONS = [];
elseif nargin > 4
        error('Too many inputs.')
end

% Set user options
opts = setQUINOPTOptions(OPTIONS);

% ----------------------------------------------------------------------- %
% Construct inequality problem
% ----------------------------------------------------------------------- %
INEQ = setInequalityModel(EXPR,BC);
[INEQ,Equalities,FLAG] = integrateByParts(INEQ);
if FLAG==1
    fprintf(['\nWARNING: Integration by parts produced an infeasible constraint.',...
            'Your problem is infeasible!\n\n']);
    varargout{3} = FLAG;
    return
end

% ----------------------------------------------------------------------- %
% Find expansions
% ----------------------------------------------------------------------- %
if ~any(INEQ.DERORD)
    % Simply set a SOS constraint, no expansion needed
    Q = 1;                                  % just set to 1 for simplicity
    S = [0.5*INEQ.F.Fb, INEQ.F.Fm; ...      % Divide Fb by 2 to "bring it inside the integral"
        zeros(size(INEQ.F.Fm))', INEQ.F.Fi];
    S = (S+S')/2;
    if ~isempty(INEQ.BC)&&~isZero(INEQ.BC)
        M = null(full(INEQ.BC));
        M = spblkdiag(M,speye(size(INEQ.F.Fi)));
        S = M'*S*M;
    end
    MatrixInequalities = {};
    slacks.t = []; slacks.pcoef = [];
    AuxVars = {};
else
    % Compute Legendre series expansion
    [Q,S,slacks,MatrixInequalities,AuxVars] = expandIntegrand(INEQ,N,opts);
end

% ----------------------------------------------------------------------- %
% Remove zero rows/cols and try to detect block-diagonal structure
% ----------------------------------------------------------------------- %
Q = makeBlkDiag(Q);
S = makeBlkDiag(S);

% ----------------------------------------------------------------------- %
% Set constraints 
% ----------------------------------------------------------------------- %
[CNSTR,DATA,FLAG] = setConstraints(INEQ.IVAR,Q,S,MatrixInequalities,slacks,...
                                  AuxVars,Equalities,FLAG,opts);

% ----------------------------------------------------------------------- %
% Set outputs 
% ----------------------------------------------------------------------- %
if nargout > 0; varargout{1} = CNSTR; end
if nargout > 1; varargout{2} = DATA; end 
if nargout > 2; varargout{3} = FLAG; end


% ----------------------------------------------------------------------- %
% Set undocumented outputs (not for general users)
% ----------------------------------------------------------------------- %


%% END CODE
end
