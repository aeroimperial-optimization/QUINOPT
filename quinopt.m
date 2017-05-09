function varargout = quinopt(varargin)

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
%     Each entry of EXPR must be a polynomial of the integration variable
%     returned by the function <a href="matlab:help('indvar')">indvar</a>, and a quadratic polynomial of the
%     dependent variables returned by the function <a href="matlab:help('depvar')">depvar</a> variables.
%     Examples can be found in the folder "examples/". The output SOL is a
%     structure with the following fields:
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
% SOL = QUINOPT(EXPR,BC,OBJ) minimizes the objective function OBJ constrained by
%     the integral inequalities specified by EXPR.
%
% SOL = QUINOPT(EXPR,BC,OBJ,OPTIONS) overrides the default options. OPTIONS is a
%     structure containing at least one of the following fields:
%
%     - YALMIP: a substructure containing the options for YALMIP, set
%               with YALMIP's command <a href="matlab:help('sdpsettings')">sdpsettings</a>.
%
%     - N: an integer specifying the number of Legendre coefficients to use in
%          the expansion of the dependent variable to obtain an
%          SDP-representable relaxation of the quadratic integral inequality.
%
%     - method: if set to 'inner' (default), QUINOPT generates an inner
%               approximation of the feasible set of the integral inequalities
%               specified by EXPR, i.e. the quadratic integral inequality is
%               strenghtened. If set to 'outer', an outer approximation is
%               constructed, i.e. the integral inequality is relaxed into a
%               weaker constraint. (NOTE: OPTIONS.method replaces the deprecated
%               option OPTIONS.rigorous).
%
%     - BCprojectorBasis: string specifying which basis to use for the
%               projection on the boundary conditions. If set to 'rref'
%               (default), use a "rational" basis. If set to 'orth', use an
%               orthonormal basis. The orthonormal basis may be preferable
%               numerically, but it may destroy sparsity of the data.
%
%     - sosdeg: the degree of the sum-of-squares polynomials used in the
%               S-procedure to localize SOS constraints from the integral
%               inequality to the integration domain. Default value: 6.
%
%     - solve: if set to 'true' (default), QUINOPT calls the solver specified by
%              the YALMIP options (or YALMIP's default solver). If set to
%              'false', QUINOPT does not call the solver, but simply sets up the
%              YALMIP problem structure. In this case, additional outputs to
%              QUINOPT are required (see below).
%
% SOL = QUINOPT(EXPR,BC,OBJ,OPTIONS,CNSTR), and 
% SOL = QUINOPT(EXPR,BC,OBJ,OPTIONS,CNSTR,PARAMETERS) minimize the objective 
%     function OBJ subjet to the integral inequalities specified by EXPR and BC,
%     and the additional constraints given by CNSTR. If CNSTR contains 
%     sum-of-square constraints, then the parameters in the polynomial 
%     expressions MUST be specified in the input vector PARAMETERS. See <a href="matlab:help('@sdpvar/sos')">sos</a> and
%     <a href="matlab:help('yalmip/modules/sos/solvesos')">solvesos</a> for more details on specifying sum-of-squares constraints with 
%     YALMIP.
%
% [SOL,CNSTR,DATA] = QUINOPT(...) also returns the YALMIP constraint object
%      CNSTR used to solve the optimization problem, and a structure DATA
%      containing all raw variables used to set up the constraints in CNSTR.
%
% Unused inputs can be left empty; for example, QUINOPT(EXPR,BC,[],CNSTR)
% determines whether the integral inequality and the constraints in CNSTR
% are feasible without optimizing an objective function.
%
% See also INDVAR, DEPVAR, QUINOPTFEASCODE, OPTIMIZE, @SDPVAR/VALUE,
%          SDPSETTINGS
%

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    05/05/2016
% Last Modified:    09/05/2017
% ----------------------------------------------------------------------- %

% Better method to clear model
if nargin==1 && ischar(varargin{1}) && strcmpi(varargin{1},'clear')
    W = evalin('caller','whos');
    for i = 1:size(W,1)
        if any(strcmp(W(i).class,{'legpoly','indvar','depvar'}))
            evalin('caller', ['clear ' W(i).name ';']);
        end
    end
    
    % Clear persistend variables
    evalin('caller','dvarlist clear');
    evalin('caller','qiimodel clear');
    return
end

% Version
if nargin==1 && ischar(varargin{1}) && strcmpi(varargin{1},'version')
    varargout{1} = '1.5';
    warning(['QUINOPT''s calling syntax has changed since version 1.5. ',...
        '<a href="matlab:help(''quinopt'')">See the help for more details</a>.'])
    return
end

% Help?
if nargin < 1;
    help('quinopt');
    return;
end

% ----------------------------------------------------------------------- %
% Get variables
% try to handle cases in which old syntax was used
EXPR = [];
BC = [];
OBJ = [];
OPTIONS = [];
CNSTR = [];
PARAMETERS = [];
N = [];
warned_user = 0;

% First three inputs unchanged
if nargin>=1; EXPR = varargin{1}; end
if nargin>=2; BC = varargin{2}; end
if nargin>=3; OBJ = varargin{3}; end

% Fourth input: a constraint (deprecated), or options structure
if nargin>=4
    if isa(varargin{4},'struct')
        OPTIONS = varargin{4};
        try
            N = OPTIONS.N;
        catch
            N = [];
        end
        
    elseif isa(varargin{4},'constraint') || isa(varargin{4},'lmi')
        if ~warned_user
            warning(['The calling syntax for quinopt() has changed since version 1.5. ',...
                'Please <a href="matlab:help(''quinopt'')">see the function help</a> for more details.'])
            warned_user = 1;
        end
        CNSTR = varargin{4};
        
    end
end

% Fifth input: an sdpvar (deprecated), or a constraint
if nargin>=5
    if isa(varargin{5},'constraint') || isa(varargin{5},'lmi')
        CNSTR = varargin{5};
        
    elseif isa(varargin{5},'sdpvar')
        if ~warned_user
            warning(['The calling syntax for quinopt() has changed since version 1.5. ',...
                'Please <a href="matlab:help(''quinopt'')">see the function help</a> for more details.'])
            warned_user = 1;
        end
        PARAMETERS = varargin{5};
        
    end
end

% Sixth input: a scalar (deprecated), or an sdpvar
if nargin>=6
    
    if isa(varargin{6},'sdpvar')
        PARAMETERS = varargin{6};
        
    elseif isnumeric(varargin{6})
        if ~warned_user
            warning(['The calling syntax for quinopt() has changed since version 1.5. ',...
                'Please <a href="matlab:help(''quinopt'')">see the function help</a> for more details.'])
            warned_user = 1;
        end
        N = varargin{6};
    end
end

% Seventh inputh (deprecated)
if nargin==7 && ( isempty(varargin{7}) || isa(varargin{7},'struct') )
    if ~warned_user
        warning(['The calling syntax for quinopt() has changed since version 1.5. ',...
            'Please <a href="matlab:help(''quinopt'')">see the function help</a> for more details.'])
    end
    OPTIONS = varargin{7};
end

% Error if more than 7 inputs
if nargin > 7; error('Too many input arguments'); end


% Compute and set outputs
[OUT1,OUT2,OUT3] = quinoptSolve(EXPR,BC,OBJ,CNSTR,PARAMETERS,N,OPTIONS);
varargout = cell(nargout,1);
if nargout > 0; varargout{1} = OUT1; end;
if nargout > 1; varargout{2} = OUT2; end;
if nargout > 2; varargout{3} = OUT3; end;
if nargout > 3; error('Too many output arguments.'); end;


% END CODE
end