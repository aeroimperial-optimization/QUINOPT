function x = indvar(a,b)

%% INDVAR.m Create independent variable of integration
%
% x = INDVAR(a,b) creates an independent variable of integration with
%       domain [a,b] used to set up an integral inequality constraint with 
%       the toolbox QUINOPT. 
%       
% NOTE: An integral inequality constraint defined with QUINOPT can only
%       have one independent variable - multivariable integrals are not
%       allowed. However, one independent variable can be used to define
%       multiple integral inequalities.
%
% See also @INDVAR/GETDOMAIN, @INDVAR/SETDOMAIN, DEPVAR, QUINOPT

% Subclass of sdpvar with additional properties:
%   - indvarID     : the unique identifier of this indvar
%   - indvarCleanup: cleanup operations to remove all dependent objects

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    25/04/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% Check domain
if ~isnumeric(a)||~isnumeric(b)||numel(a)~=1||numel(b)~=1||~isreal(a+b)
    error('Inputs a and b must be real numeric scalars.')
elseif any(isinf([a,b]))
    error('The domain for the independent variable must be bounded.')
end

% Create variable
y = sdpvar(1,1); 

% Set domain in model
qiimodel('indvar',y,a,b);

% Set class as a subclass of sdpvar y
x.indvarCleanup = onCleanup(@()delete(depends(y)));
x = class(x,'indvar',y);

end

%% NESTED DESTRUCTOR
function delete(IVARID)
    % Find and remove dependent variables from internal model
    mod = qiimodel('query');
    if ~isempty(mod.DEPVARMODEL)
        dvarpos = find(depends(mod.DEPVARMODEL.IVAR)==IVARID);
        b = getbase(mod.DEPVARMODEL.IVAR);
        dvarpos = find(b(:,dvarpos+1))';
        if ~isempty(dvarpos)
            for i = 1:length(dvarpos)
                qiimodel('cleardepvar',mod.DEPVARMODEL.DVARID(dvarpos(i)));
            end
        end
    end

    % Clear independent variable from internal model (if not already cleared?)
    if ~isempty(mod.INDVARMODEL)
        qiimodel('clearindvar',IVARID);
    end

end
