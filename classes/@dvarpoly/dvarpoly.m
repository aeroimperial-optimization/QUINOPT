function dvar = dvarpoly(varargin)

%% DVARPOLY.m polynomial in dependent variables
%
% Define a class to represent polynomials in a convenient way, compatible
% with YALMIP and the legendre polynomials of class <a href="matlab:help('legpoly')">legpoly</a>.
%
% Properties of each class instance:
%        -coeff       List of coefficients of each monomial
%        -ivars       List of independent variables. [] means no variable 
%                     (a constant polynomial).
%        -monom       Matrix describing monomials. Each row corresponds to
%                     a monomial. Entries of each row represent powers of
%                     corresponding independent variable listed in ivars.
%
% See also DEPVAR

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    16/04/2015
% Last Modified:    16/04/2016
% ----------------------------------------------------------------------- %

% CODE
superiorto('double');
superiorto('sdpvar');
superiorto('indvar');
superiorto('legpoly');

% Short syntax for scalar variables
if ischar(varargin{1})
    % Construct variables with given name
    for k = 1:nargin
        if ~isvarname(varargin{k})
            error('Invalid variable name.')
        else
            assignin('caller',varargin{k},eval('dvarpoly(1,1);'));
        end
    end
    return
end


% Proper syntax for vectors  & converters

switch nargin
    
    % ------------------------------------------------------------------- %
    case 0
        % New variable
        dvar = dvarpoly(1,1);
        return
        
        
        % ------------------------------------------------------------------- %
    case 1
        % Class converters
        
        if isstruct(varargin{1})
            % Convert structure input to depvarpoly
            try
                [m,n] = size(varargin{1});
                [dvar(1:m,1:n).coeff] = deal(varargin{1}.coeff);
                % [dvar(1:m,1:n).ivars]= deal(varargin{1}.ivars(:));
                [dvar(1:m,1:n).ivars]= deal(varargin{1}.ivars);
                [dvar(1:m,1:n).monom] = deal(varargin{1}.monom);
                dvar = class(dvar,'dvarpoly');
                return
                
            catch
                error('Input structure does not have the required fields.')
                
            end
            
            
        elseif isnumeric(varargin{1}) || isa(varargin{1},'sdpvar') || isa(varargin{1},'legpolyobj')
            % Converter a double, sdpvar, legpolyobj
            dvar.coeff = varargin{1};
            dvar.ivars = [];
            dvar.monom = 0;
            dvar = class(dvar,'dvarpoly');
            return
            
        end
        
        % ------------------------------------------------------------------- %
    case 2
        % Get size of output
        m = varargin{1};
        n = varargin{2};
        
        % Update list of existing variables
        VLIST = dvarlist;
        newvars = VLIST(end)+1:VLIST(end)+n*m;
        dvarlist(newvars);
        newvars = num2cell(newvars);
        
        % Set output
        [dvar(1:m,1:n).coeff] = deal(1);
        [dvar(1:m,1:n).ivars] = deal(newvars{:});
        [dvar(1:m,1:n).monom] = deal(1);
        dvar = class(dvar,'dvarpoly');
        return
        
    otherwise
        error('Too many input arguments.')
        
end

% ----------------------------------------------------------------------- %
% END CODE
end
