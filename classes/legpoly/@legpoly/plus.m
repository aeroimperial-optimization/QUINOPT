function s = plus(X,Y)

% plus.m
%
% c = plus(X,Y)

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

% CODE

% ----------------------------------------------------------------------- %
% LEGPOLY + NUMERIC
% ----------------------------------------------------------------------- %
if isa(X,'legpoly') && isnumeric(Y)
    
    [m,n] = size(X);
    
    if size(X)==size(Y)
        % Define element-wise addition of a constant
        
        [s(1:m,1:n).coef] = deal(X.coef);
        newcoeff = arrayfun(@(x,y)[x.coef(1)+y;x.coef(2:end)],s,Y,'uniformoutput',0);
        [s(1:m,1:n).coef] = deal(newcoeff{:});
        [s(1:m,1:n).ivar] = deal(X.ivar); 
        [s(1:m,1:n).domn] = deal(X.domn); 
        s = legpoly(s);
        return
        
    elseif numel(Y)==1
        % Inflate Y to do element-wise addition of a constant
        Y = repmat(Y,m,n);
        s = X+Y;
        
    else
        error('Dimensions of addends mismatch.')
        
    end
    
elseif isnumeric(X) && isa(Y,'legpoly')
    s = Y+X;
    return
    
% ----------------------------------------------------------------------- %
% LEGPOLY + LEGPOLY
% ----------------------------------------------------------------------- %
elseif isa(X,'legpoly') && isa(Y,'legpoly')
    % Addition between 2 legpolys
    [m,n] = size(X);
    ivarX = getivar(X);
    ivarY = getivar(Y);
    domX = getDomain(X);
    domY = getDomain(Y);
    degX = degree(X);
    degY = degree(Y);
    
    if depends(ivarX)~=depends(ivarY)
        % if Y has degree 0, ivarY is the scalar 0 rather than an sdpvar.
        % However, depends(0)=[] and A~=[] returns [], which is interpreted
        % as false. This is good, since can add a constant to any legpoly!
        error(['Legendre polynomials have different independent variables. '...
            'Only operations between polynomials with the same independent variable are supported.'])
    
    elseif any(domX-domY) && any(domX) && any(domY)
        error('Cannot add Legendre polynomials defined over different intervals.')
        
    elseif ~any(size(X)-size(Y)) %size(X)==size(Y)
        % Define element-wise addition of a legpolys
        
        % Find cumulative independent variables
        if degX ~=0
            index = find(ivarX);
            ivar = num2cell(ivarX(index(1)).*any(ivarX+ivarY));
        elseif degY ~=0
            index = find(ivarY);
            ivar = num2cell(ivarY(index(1)).*any(ivarX+ivarY));
        else
            [ivar{1:m,1:n}] = deal(0);
        end
        
        % Structure with data from X
        degx = num2cell(degree(X,[],'all'));
        [xstr(1:m,1:n).deg] = deal(degx{:});
        [xstr(1:m,1:n).coef] = deal(X.coef);
        
        % Structure with data from Y
        degy = num2cell(degree(Y,[],'all'));
        [ystr(1:m,1:n).deg] = deal(degy{:});
        [ystr(1:m,1:n).coef] = deal(Y.coef);
        
        % Find new coefficients
        newdeg = arrayfun(@(x,y)max(x.deg,y.deg),xstr,ystr,'uniformoutput',0);
        [s(1:m,1:n).deg] = deal(newdeg{:});
        newcoeff = arrayfun(@(x,y,z)[x.coef; zeros(z.deg-x.deg,1)] ...
                                  + [y.coef; zeros(z.deg-y.deg,1)],xstr,ystr,s,'uniformoutput',0);
        newcoeff = cellfun(@(x)removeTrailingZeros(x),newcoeff,'uniformoutput',0);
        [s(1:m,1:n).coef] = deal(newcoeff{:});
        
        % Update independent variable
        ivar = cellfun(@(x,y)(length(x)>1)*y,newcoeff,ivar,'uniformoutput',0);
        [s(1:m,1:n).ivar] = deal(ivar{:});
        
        % Update domain
        domn = cellfun(@(x,y)any(x)*domX,ivar,'uniformoutput',0);
        [s(1:m,1:n).domn] = deal(domn{:});
        
        % Set output
        s = legpoly(s);
        return
    
    elseif numel(X)==1
        X = repmat(X,size(Y,1),size(Y,2));
        s = X+Y;
        return
   
    elseif numel(Y)==1
        Y = repmat(Y,size(X,1),size(X,2));
        s = X+Y;
        return
    
    else
        error('Dimensions of addends mismatch.')
        
    end
    
% ----------------------------------------------------------------------- %
% LEGPOLY + SDPVAR
% ----------------------------------------------------------------------- %
elseif isa(X,'legpoly') && isa(Y,'sdpvar')
    
    % Convert Y to legpoly & sum
    ivarX = getivar(X);
    idx = find(any(ivarX),1,'first');
    domn = getDomain(X);
    if ~isempty(idx)
        ivarX = X(idx).ivar;
        Y = legpoly(domn,ivarX,Y);
    else
    	% X is constant polynomial - better to transform X!
        X = coefficients(X);
    end
    s = X+Y;
	s = legpoly(s);
    return
    
elseif isa(X,'sdpvar') && isa(Y,'legpoly')
    s = Y+X;
    return
    
    % ----------------------------------------------------------------------- %
    % END SCRIPT
end

