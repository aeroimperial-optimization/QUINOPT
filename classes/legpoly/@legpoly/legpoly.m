function [P,C] = legpoly(varargin)

%% LEGPOLY Polynomial in Legendre basis (part of the QuadIntIneq toolbox)
%
% P = LEGPOLY(x,DEG) creates a polynomial P in the independent variable x
%       of degree DEG, expressed in Legendre basis. That is, P(x) is
%       expressed as
%
%       P(x) = C(1)*L_0[z(x)] + C(2)*L_1[z(x)] + ... + C(DEG+1)*L_DEG[z(x)]
%
%       where L_n(z) is the Legendre polynomial of degree n. Since Legendre
%       polynomials are defined over the standard domain [-1,1], the
%       original independent variable x with domain [a,b] is rescaled to
%
%           z(x) = (2*x-b-a)/(b-a)
%
%       The input x must be a valid independent variable (class indvar),
%       and DEG should be a non-negative integer. The coefficients of the
%       polynomial are YALMIP variables (sdpvar objects) and can be
%       recovered with the command "C = coefficients(P)". Finally, P can be 
%       displayed symbolically in the standard monomial basis using the 
%       command "sdisplay(P)".
%
% [P,C] = LEGPOLY(x,DEG) also returns the Legendre coefficients of the
%       polynomials in the vector C. These are YALMIP variables (sdpvar
%       objects). The coefficients in C are listed in order of increasing
%       degree of the corresponding Legendre polynomial (see above).
%
% P = LEGPOLY(x,DEG,COEF) creates a polynomial P expressed in Legendre
%       basis whose coefficients are specified by COEF. COEF can be a
%       numeric/sdpvar vector, or an M-by-N cell array whose entries are
%       numeric/sdpvar vectors. When COEF is a cell array, an M-by-N matrix
%       of Legendre polynomials is created such that the coefficients of
%       the entry P(i,j) are given by COEF{i,j}.
%
% [P,C] = LEGPOLY(x,DEG,ROWS,COLS), creates a ROWS-by-COLS matrix of
%       Legendre polynomials of degree DEG. The output C is optional.
%
% EXAMPLE.
%
% >> x = indvar(0,1);       % define independent variable with domain [0,2].
% >> P = legpoly(x,2);      % define polynomial in Legendre basis of degree 2.
% >> C = coefficients(P);   % recover the variable coefficients of P.
% >> sdisplay(P);           % print P to screen (in monomial basis).
%
% See also INDVAR, @LEGPOLY/COEFFICIENTS, @LEGPOLY/SDISPLAY


% ----------------------------------------------------------------------- %
% ADVANCED INTERNAL USES FOR CLASS CONVERSION
%
% P = LEGPOLY(N) converts the numeric/sdpvar scalar N to a constant Legendre
%       polynomial and is equivalent to P = LEGPOLY(x,0,N).
%
% P = LEGPOLY(x,P) converts an N-by-M matrix of sdpvar polynomials P to an
%       N-by-M matrix of polynomials in the Legendre basis with independent
%       variable x.
%
% P = LEGPOLY(dom,x,P) converts the polynomial P to a legpoly when the
%       independent variable is an sdpvar but not an independent variable 
%       (class indvar). In this case, the domain of the independent variable must be
%       specified by an additional input dom.
%
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% Properties of a legpoly P
% P.coef:   the Legendre coefficients (a vector)
% P.ivar:   the independent variable (an sdpvar or 0 for constants)
% P.domn:   the domain of the independent variable (a vector [a b] or [0 0] for constants)

superiorto('double');
superiorto('sdpvar');
superiorto('indvar');

switch nargin
    
    case 0
        % ZERO LEGPOLY
        P = struct('coef',0,'ivar',0,'domn',[0,0]);
        P = class(P,'legpoly');
        return
    
    % ------------------------------------------------------------------- %
    case 1
    % ------------------------------------------------------------------- %
    % CLASS CONVERTERS (internal use)
        
        x = varargin{1};
        
        if isstruct(x) && nargout<=1
            % Assign class legpoly to input structure x
            [m,n] = size(x);
            try
                % Set coefficients removing trailing zeros
                C = arrayfun(@(z)removeTrailingZeros(z.coef),x,'uniformoutput',0);
                [P(1:m,1:n).coef] = deal(C{:});
                
                % Set independent variable checking for degree 0
                ivar = arrayfun(@(z)(length(z.coef)>1)*z.ivar,x,'uniformoutput',0);
                [P(1:m,1:n).ivar] = deal(ivar{:});
                
                % Set domain of independent variable checking for degree 0
                domn = arrayfun(@(z)(length(z.coef)>1)*z.domn,x,'uniformoutput',0);
                [P(1:m,1:n).domn] = deal(domn{:});
                
                % Assign class
                P = class(P,'legpoly');
                return
                
            catch
                error('Input structure does not have the required fields.')
            end
            
        elseif (isnumeric(x)||isa(x,'sdpvar')) && nargout<=1
            % Convert vector of double/sdpvar to vector of legpolys of degree 0
            [m,n] = size(x);
            
            % Set coefficients
            C = x(:);
            coef = num2cell(C);
            [P(1:m,1:n).coef] = deal(coef{:});
            
            % Set independent variable to 0
            [P(1:m,1:n).ivar] = deal(0);
            
            % Set independent variable domain to [0,0] (identifies a constant)
            [P(1:m,1:n).domn] = deal([0,0]);
            
            % Assign class
            P = class(P,'legpoly');
            return
            
        else
            error('At least two input arguments are required. Type "help legpoly" for help.')
        end
        
    % ------------------------------------------------------------------- %
    case 2
    % ------------------------------------------------------------------- %
    % CONSTRUCT 1-by-1 POLYNOMIAL OF DEGREE varargin{2}
        
        % Check inputs
        x = varargin{1};
        if  ~isa(x,'indvar')||~isValidindvar(x)
            error('Input x must be a valid independent variable. Type "help indvar" for more information.')
        end
        
        % Construct object
        DOMAIN = getDomain(x);

        if isnumeric(varargin{2}) && isscalar(varargin{2})
            % Create 1-by-1 polynomial of degree varargin{2}
            if varargin{2}<0
                error('Cannot have negative polynomial degree.')
            elseif varargin{2}>0
                C = sdpvar(varargin{2}+1,1);
                P.coef = C;
                P.ivar = sdpvar(x);
                P.domn = DOMAIN;
            else
                C = sdpvar(varargin{2}+1,1);
                P.coef = C;
                P.ivar = sdpvar(x);
                P.domn = DOMAIN;
            end
            
            % Assign class
            P = class(P,'legpoly');
            return

        elseif degree(varargin{2},x)~=0
            % vector of sdpvar polynomials to convert
            % VERY UGLY LOOP!
            S = varargin{2};
            [m,n] = size(S);
            C = cell(m,n);
            for i = 1:m
                for j = 1:n
                    deg = degree(S(i,j),x);      % sdpvar/degree
                    if deg == 0
                        C{i,j} = S(i,j);
                        P(i,j).coef = S(i,j);
                        P(i,j).ivar = 0;
                        P(i,j).domn = [0,0];
                    else
                        % Replace variable
                        z = ( (DOMAIN(2)-DOMAIN(1))*x + DOMAIN(2)+DOMAIN(1) )/2;
                        q = replace(S(i,j),x,z);
                        LegCoef = legBasisCoef(q,x);
                        C{i,j} = LegCoef(:);
                        P(i,j).coef = LegCoef(:);
                        P(i,j).ivar = x;
                        P(i,j).domn = DOMAIN;
                    end
                end
            end
            
            % Assign class
            if m*n==1; C = C{1}; end
            P = class(P,'legpoly');
            return
        
        else
            error('Input DEG must be a non-negative integer.')

        end
        
    % ------------------------------------------------------------------- %
    case 3
    % ------------------------------------------------------------------- %
   
        if isa(varargin{1},'sdpvar')
            % CONSTRUCT POLYNOMIAL WITH GIVEN COEFFICIENTS (IGNORE DEGREE INPUT BUT
            % DISPLAY WARNING IF MISMATCHING)

            % Check inputs
            x = varargin{1};
            if  ~isa(x,'indvar')||~isValidindvar(x)
                error('Input x must be a valid independent variable. Type "help indvar" for more information.')
            elseif ~isnumeric(varargin{2}) || ~isscalar(varargin{2})
                error('Input DEG must be a suitable positive integer.')
            end
            DOMAIN = getDomain(x);
            x = sdpvar(x);
            
            if ~iscell(varargin{3}) && isvector(varargin{3})
                % input is a vector of Legendre coefficients to be used to
                % construct a legpoly
                if length(varargin{3})~=varargin{2}+1
                    warning(['Specified degree does not match the number of coefficients provided. ',...
                        'I''ll ignore the input degree and use the coefficients instead.'])
                end
                
                C = removeTrailingZeros(varargin{3}(:));
                if ~isZero(C)
                    P.coef = C;
                    if numel(C)>1
                        P.ivar = x;
                        P.domn = DOMAIN;
                    else
                        % Degree 0
                        P.ivar = 0;
                        P.domn = [0,0];
                    end
                else
                    % The 0 polynomial
                    P.coef = 0;
                    P.ivar = 0;
                    P.domn = [0,0];
                end
                
                % Assign class
                P = class(P,'legpoly');
                return
                
            elseif iscell(varargin{3})
                % create m-by-n polynomial with coefficients of
                % entry (i,j) are specified by varargin{3}(i,j).
                
                C = varargin{3};
                [m,n] = size(C);
                
                % Check that all entries in cell are vectors
                notvec = ~cellfun(@isvector,C(:));
                if nnz(notvec)~=0
                    error('All entries of the cell array of coefficients must contain a non-empty vector.')
                end
                
                % Compute coefficients and find independent variable & domain
                C = cellfun(@removeTrailingZeros,C(:),'uniformoutput',0);
                ivar = cellfun(@(y)(length(y)>1)*x,C,'uniformoutput',0);
                domn = cellfun(@(y)(length(y)>1)*DOMAIN,C,'uniformoutput',0);
                
                % Assign fields
                [P(1:m,1:n).coef] = deal(C{:});
                [P(1:m,1:n).ivar] = deal(ivar{:});
                [P(1:m,1:n).domn] = deal(domn{:});
                
                % Assign class
                P = class(P,'legpoly');
                return
                
            else
                error('Third input to "legpoly" must be a scalar, a vector or a cell array of vectors.')
                
            end
        else
            % CONVERT POLYNOMIAL OF SDPVAR ONLY - NO INDVAR (INTERNAL USE)
            DOMAIN = varargin{1};
            x = varargin{2};
            S = varargin{3};
            [m,n] = size(S);
            C = cell(m,n);
            for i = 1:m
                for j = 1:n
                    deg = degree(S(i,j),x);      % sdpvar/degree
                    if deg == 0
                        C{i,j} = S(i,j);
                        P(i,j).coef = S(i,j);
                        P(i,j).ivar = 0;
                        P(i,j).domn = [0,0];
                    else
                        % Replace variable
                        z = ( (DOMAIN(2)-DOMAIN(1))*x + DOMAIN(2)+DOMAIN(1) )/2;
                        q = replace(S(i,j),x,z);
                        LegCoef = legBasisCoef(q,x);
                        C{i,j} = LegCoef(:);
                        P(i,j).coef = LegCoef(:);
                        P(i,j).ivar = x;
                        P(i,j).domn = DOMAIN;
                    end
                end
            end
            
            % Assign class
            if m*n==1; C = C{1}; end
            P = class(P,'legpoly');
            return

        end
        
    % ------------------------------------------------------------------- %
    case 4
    % ------------------------------------------------------------------- %
    % CREATE ROWS-by-COLS POLYNOMIAL OF DEGREE varargin{3}
        
        % Check if first inputs is a suitable sdpvar
        x = varargin{1};
        if  ~isa(x,'indvar')||~isValidindvar(x)
            error('Input x must be a valid independent variable. Type "help indvar" for more information.')
        end
        DOMAIN = getDomain(x);
        x = sdpvar(x);
        deg = varargin{2};
        m = varargin{3};
        n = varargin{4};
        if  ~isnumeric(deg) || numel(deg)~=1 || rem(deg,1)~=0
            error('Specified degree must be a scalar integer.')
        elseif  ~isnumeric(m) || ~isscalar(m) || ~isnumeric(n) || ~isscalar(n) || m*n<=0
            error('Inputs ROWS and COLS must be positive integers.')
        end
        
        
        % Create coefficients
        C = cell(m*n,1);
        for i=1:numel(C)
            C{i} = sdpvar(deg+1,1);
        end
        
        % Find independent variable & domain based on degree
        ivar = cellfun(@(y)(length(y)>1)*x,C,'uniformoutput',0);
        domn = cellfun(@(y)(length(y)>1)*DOMAIN,C,'uniformoutput',0);
        
        % Assign fields
        [P(1:m,1:n).coef] = deal(C{:});
        [P(1:m,1:n).ivar] = deal(ivar{:});
        [P(1:m,1:n).domn] = deal(domn{:});
        
        % Assign class
        C = reshape(C,m,n);
        P = class(P,'legpoly');
        return
        
end


% END FUNCTION - CLASS CONSTRUCTOR
end