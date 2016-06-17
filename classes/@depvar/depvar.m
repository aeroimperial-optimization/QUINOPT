function varargout = depvar(x)

%% DEPVAR.m Define dependent variables
%
% U = DEPVAR(x) sets up a symbolic variable modelling a generic function U(x), 
%       where x is a valid independent variable with domain [a,b] created 
%       with the command <a href="matlab:help('indvar')">indvar</a>. U behaves like a function handle and is
%       used with the syntax 
%
%       U(POINT) 
%
%       where POINT is either the independent variable x, the lower extremum
%       a of the domain of U, or the upper extremum b of the domain of U.
%       Moreover, derivatives of U can be created using the syntax
%
%       U(POINT,DERIVATIVE)
%
%       where POINT is x, a or b and DERIVATIVE is the desired derivative.
%       Note that U(POINT,0) is equivalent to U(POINT). See below for
%       some examples.
%     
% [U1,U2,...Uq] = DEPVAR(x)  sets up multiple dependent variables. Each
%       dependent variable depends on the independent variable x. 
%
% EXAMPLES.
%
% >> x = indvar(0,1);       % setup independent variable with domain [0,1]
% >> [u,v] = DEPVAR(x);     % setup dependent variables u, v.
% >> u(x)                   % returns the symbolic variable used to model u(x)
% >> v(1)                   % returns the symbolic variable used to model the 
%                           % boundary value v(1)
% >> v(x,2)                 % returns the symbolic variable used to model the 
%                           % second derivative v''(x)
% >> u(0,1)                 % returns the symbolic variable used to model the 
%                           % boundary value of first derivative u'(0)
% >> v(0.5)                 % returns an error since POINT is different from
%                           % the independent variable x or the boundary 
%                           % points 0 and 1
% >> P = u(x,1)^2-v(0)^2    % constructs a polynomial P of dependent variables 
%                           % and boundary values.
%
% NOTE: dependent variables created using DEPVAR become invalid when the
%       independent variable x on which they depend is cleared from the
%       workspace, but are not deleted.
%
% See also INDVAR, QUINOPT, SETQUADINTINEQ, QIIMODEL, DVARPOLY

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    25/04/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

%% MAIN CLASS DEFINITION
superiorto('sdpvar')

% Set default arguments
if nargin < 1
    error('Not enough input arguments.');
elseif nargin > 1
    error('Too many input arguments.');
end

% Initalise output
ndvars = max(nargout,1);
varargout = cell(ndvars,1);

% Call constructors
for i=1:ndvars
    
    % Check inputs
    if ~isa(x,'indvar')||~isValidindvar(x)
        error('Input x must be a valid scalar independent variable.')
    end
    
    % Set dependent variables
    dvarID = qiimodel('depvar',x,0);
    varargout{i}.id = dvarID;
    varargout{i}.h = @(POINT,DER) dependentVariableHandle(dvarID,POINT,DER);
    varargout{i}.depvarCleanup = onCleanup(@()delete(dvarID));
    varargout{i} = class(varargout{i},'depvar');
    
end

end


%% NESTED FUNCTION TO BUILD ANONYMOUS FUNCTION HANDLE
% The following cannot be saved as a separate file otherwise the handle
% would not be anonymous and would have problems in calling with different
% input types.

function u = dependentVariableHandle(DVARID,POINT,DER)
    %% DEPENDENTVARIABLEHANDLE.m Return dependent variables from internal model
    %
    % u = DEPENDENTVARIABLEHANDLE(DVARID) returns the internal variable used by
    %       QUINOPT to model the basic dependent variable with identifier
    %       DVARID. That is, the zeroth derivative is returned.
    %
    % u = DEPENDENTVARIABLEHANDLE(DVARID,POINT) returns the internal variable
    %       used by QUINOPT to model the basic dependent variable with
    %       identifier DVARID evaluated at the point POINT. POINT can be:
    %
    %       - the independent variable on which the dependent variable depends
    %       - one of the boundary values of the independent variable.
    %
    %       The allowed boundary values and the correct independent variables are
    %       determined from the internal model. If POINT differs from the
    %       correct independent variable or the allowed boundary values, an
    %       error is thrown.
    %
    % u = DEPENDENTVARIABLEHANDLE(DVARID,POINT,DER) returns the internal
    %       variable used by QUINOPT to model the derivative of order DER
    %       of the dependent function with identifier DVARIND, evaluated at
    %       point POINT. If variables modelling the requested derivatives do
    %       not already exist, new internal variables to the internal model.
    %
    % NOTE: This function cannot be used on its own - it must be defined as
    %
    % See also DEPVAR, INDVAR, QIIMODEL

    % Initial check
    if nargin==0
        error('Not enough input argument: dependent variable identifier must be specified.')
    end

    % Extract model for this variable
    mod = qiimodel('query');
    if isempty(mod.DEPVARMODEL)
        error('No dependent variable exists. Possible cause: you cleared all independent variable.')
    end
    varInd = find(mod.DEPVARMODEL.DVARID==DVARID);
    if isempty(varInd)
        error('Invalid dependent variable. Possible cause: you cleared the corresponding independent variable.')
    end
    IVAR = mod.DEPVARMODEL.IVAR(varInd);
    IVARind = find(mod.INDVARMODEL.IVARID==depends(IVAR),1);
    MAXDER = mod.DEPVARMODEL.MAXDER(varInd);

    % Set optional inputs to default if not given
    if nargin==1
        POINT = IVAR;
        DER = 0;
    elseif nargin==2
        DER = 0;
    end

    % Only one point and one derivative?
    if numel(POINT)~=1 || numel(DER)~=1
        error('You can only ask for one derivative at one point.')
    end

    % Do we need to add variables & update model?
    if DER > MAXDER
        qiimodel('adddepvarderivative',DVARID,DER);
        mod = qiimodel('query');
    end

    % Return output
    firstDVARind = mod.DEPVARMODEL.DVARPART(varInd)+1;
    firstBVALind = mod.DEPVARMODEL.BVALPART(varInd)+1;
    DOMAIN = mod.INDVARMODEL.DOMAIN(IVARind,:);
    if isa(POINT,'sdpvar')
        if depends(POINT)==depends(IVAR)
            u = mod.DEPVARMODEL.DVAR(firstDVARind+DER);
        else
            error('The requested dependent variable does not depend on the specified independent variable.')
        end

    elseif POINT == DOMAIN(1)
        u = mod.DEPVARMODEL.BVAL(firstBVALind+2*DER);

    elseif POINT == DOMAIN(2)
        u = mod.DEPVARMODEL.BVAL(firstBVALind+2*DER+1);

    else
        error('Input POINT must be the independent variable or the domain endpoints (%.6g or %.6g).',...
            DOMAIN(1),DOMAIN(2))
    end

end

%% NESTED CLASS DESTRUCTOR
% Clean up internal model before clearing dependent variable from workspace.
function delete(id)
    mod = qiimodel('query');
    if ~isempty(mod.DEPVARMODEL) 
        qiimodel('cleardepvar',id);
    end
end

