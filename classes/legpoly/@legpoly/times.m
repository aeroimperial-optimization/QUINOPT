function s = times(X,Y)

%% times.m
%
% s = times(X,Y)

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
%   Description:    Multiply two legpolys.
% ----------------------------------------------------------------------- %

%% CODE

% CHECK ON INPUTS
[m,n] = size(X);
[my,ny] = size(Y);
if numel(X)~=1 && numel(Y)~=1 && (m~=my || n~=ny)
    error('Size mismatch.')
end

% Any zero inputs for some reason? Need to maintain class :legpoly input =>
% legpoly output
if isZero(X) && ~isnumeric(X)
    s = Y.*0;
    if ~isa(s,'legpoly')
        s = legpoly(s);
    end
    return
    
elseif isZero(Y) && ~isnumeric(Y)
    s = X.*0;
    if ~isa(s,'legpoly')
        s = legpoly(s);
    end
    return
    
end

% ----------------------------------------------------------------------- %
% LEGPOLY x NUMERIC
% ----------------------------------------------------------------------- %
if isa(X,'legpoly') && isnumeric(Y)
    
    % Make full if sparse
    if issparse(Y) 
        Y = full(Y); 
    end
    
    if size(X)==size(Y)
        % Define element-wise multiplication
        
        % new coefficients
        [s(1:m,1:n).coef] = deal(X.coef);
        newcoeff = arrayfun(@(x,y)y*x.coef,s,Y,'uniformoutput',0);
        newcoeff = cellfun(@(x)removeTrailingZeros(x),newcoeff,'uniformoutput',0);
        [s(1:m,1:n).coef] = deal(newcoeff{:});
        
        % Update independent variable
        ivar = num2cell(getivar(X));
        ivar = cellfun(@(x,y)(length(x)>1)*y,newcoeff,ivar,'uniformoutput',0);
        [s(1:m,1:n).ivar] = deal(ivar{:});
        
        % Update domain
        domn = getDomain(X);
        domn = cellfun(@(x,y)(length(x)>1)*domn,newcoeff,'uniformoutput',0);
        [s(1:m,1:n).domn] = deal(domn{:});
        
        % Set output
        s = legpoly(s);
        return
        
    elseif numel(Y)==1
        % Inflate Y to do element-wise multiplication
        Y = repmat(Y,m,n);
        s = X.*Y;
        return
        
    elseif numel(X)==1
        % Inflate X to do element-wise multiplication
        X = repmat(X,my,ny);
        s = X.*Y;
        return
        
    else
        error('Dimensions of factors mismatch.')
    end
    
elseif isnumeric(X) && isa(Y,'legpoly')
    s = Y.*X;
    return
    
% ----------------------------------------------------------------------- %
% LEGPOLY x SDPVAR
% ----------------------------------------------------------------------- %
elseif isa(X,'legpoly') && isa(Y,'sdpvar')
    
    % Get independent variable and domain of legpoly
    ivarX = getivar(X);
    idx = find(any(ivarX),1,'first');
    domn = getDomain(X);
    
    if ~isempty(idx)
        ivarX = X(idx).ivar;
        Ydeg = degree(Y,ivarX);
        if Ydeg~=0
            % Convert Y to legpoly and multiply
            Y = legpoly(domn,ivarX,Y);
            s = X.*Y;
            return
            
        elseif size(X)==size(Y)
            % Define element-wise multiplication by an sdpvar scalar
            
            yfact = num2cell(Y(:));
            [s(1:m,1:n).coef] = deal(X.coef);
            [z(1:m,1:n).factor] = deal(yfact{:});
            newcoeff = arrayfun(@(x,y)y.factor*x.coef,s,z,'uniformoutput',0);
            newcoeff = cellfun(@(x)removeTrailingZeros(x),newcoeff,'uniformoutput',0);
            [s(1:m,1:n).coef] = deal(newcoeff{:});
        
            % Update independent variable
            ivar = num2cell(getivar(X));
            ivar = cellfun(@(x,y)(length(x)>1)*y,newcoeff,ivar,'uniformoutput',0);
            [s(1:m,1:n).ivar] = deal(ivar{:});
        
            % Update domain
            domn = getDomain(X);
            domn = cellfun(@(x,y)(length(x)>1)*domn,newcoeff,'uniformoutput',0);
            [s(1:m,1:n).domn] = deal(domn{:});
            
            % Set output
            s = legpoly(s);
            return
            
        elseif numel(Y)==1
            % Inflate Y to do element-wise multiplication
            Y = repmat(Y,m,n);
            s = X.*Y;
            return
            
        else
            error('Dimensions of factors mismatch. Only element-wise multiplications with legpolys are allowed.')
        end
        
    else
        % X is constant polynomial - better to transform X!
        X = coefficients(X);
        s = X.*Y;
        return
    end
    
 
elseif isa(X,'sdpvar') && isa(Y,'legpoly')
    s = Y.*X;
    return
    
    
% ----------------------------------------------------------------------- %
% LEGPOLY x LEGPOLY
% ----------------------------------------------------------------------- %
elseif isa(X,'legpoly') && isa (Y,'legpoly')
    % rather ugly code
    
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
        
    elseif any(domX-domY)
        error('Cannot multiply Legendre polynomials defined over different intervals.')
        
    elseif ~any(size(X)-size(Y)) %size(X)==size(Y)
        % Define element-wise product of a legpolys
        
        % Find cumulative independent variables
        % "find" for sdpvar cannot return only first index
        if degX ~=0
            index = find(ivarX);
            % ivar = num2cell(ivarX(index(1)).*any(ivarX+ivarY));
            ivar = ivarX(index(1));
        elseif degY ~=0
            index = find(ivarY);
            % ivar = num2cell(ivarY(index(1)).*any(ivarX+ivarY));
            ivar = ivarY(index(1));
        else
            % ivar = {0};
            ivar = 0;
        end
        
        % UGLY CODE
        s(m,n) = struct('ivar',[],'coef',[],'domn',[]);
        for i = 1:m
            for j = 1:n
                
                degx = degree(X(i,j));
                degy = degree(Y(i,j));
                newdeg = degx+degy;
                if newdeg==0
                    % Easy case - constants!
                    s(i,j).ivar = 0;
                    s(i,j).coef = X(i,j).coef.*Y(i,j).coef;
                    s(i,j).domn = [0,0];
                    
                elseif degx==0 || degy==0
                    % Easy case - multiply polynomial by a constants!
                    % X.coef or Y.coef is a scalar
                    s(i,j).ivar = ivar;
                    s(i,j).coef = X(i,j).coef.*Y(i,j).coef;
                    s(i,j).domn = domX;
                    
                else
                    % OLD: convert to sdpvar, multiply and convert back
%                     px = sdpvar(X(i,j));
%                     py = sdpvar(Y(i,j));
%                     q = legpoly(domX,ivar,px*py);
                    % NEW, 17/06/2016: multiplication in legendre basis
                    s(i,j).ivar = ivar;
                    s(i,j).coef = legProd(X(i,j).coef,Y(i,j).coef);
                    s(i,j).domn = domX;
                    
                end
                
            end
        end
        
        % Assign class
        s = legpoly(s);
        return
        
    elseif numel(X)==1
        % Product by a scalar
        
        % Find cumulative independent variables
        % "find" for sdpvar cannot return only first index
        if degX ~=0
            index = find(ivarX);
            % ivar = num2cell(ivarX(index(1)).*any(ivarX+ivarY));
            ivar = ivarX(index(1));
        elseif degY ~=0
            index = find(ivarY);
            % ivar = num2cell(ivarY(index(1)).*any(ivarX+ivarY));
            ivar = ivarY(index(1));
        else
            % ivar = {0};
            ivar = 0;
        end
        
        % UGLY CODE
        s(my,ny) = struct('ivar',[],'coef',[],'domn',[]);
        for i = 1:my
            for j = 1:ny
                degy = degree(Y(i,j));
                newdeg = degX+degy;
                if newdeg==0
                    % Easy case - constants!
                    s(i,j).ivar = 0;
                    s(i,j).coef = X.coef.*Y(i,j).coef;
                    s(i,j).domn = [0,0];
                    
                elseif degX==0 || degy==0
                    % Easy case - multiply polynomial by a constants!
                    % X.coef or Y.coef is a scalar
                    s(i,j).ivar = ivar;
                    s(i,j).coef = X.coef.*Y(i,j).coef;
                    s(i,j).domn = domX;
                    
                else
                    % OLD: convert to sdpvar, multiply and convert back
%                     px = sdpvar(X(i,j));
%                     py = sdpvar(Y(i,j));
%                     q = legpoly(domX,ivar,px*py);
                    % NEW, 17/06/2016: multiplication in legendre basis
                    s(i,j).ivar = ivar;
                    s(i,j).coef = legProd(X.coef,Y(i,j).coef);
                    s(i,j).domn = domX;
                    
                end
                
            end
        end
        
        % Assign class
        s = legpoly(s);
        return
        
        
        
    elseif numel(Y)==1
        % Reuse previous case (minimum overhead)
        s = Y.*X;
        return
        
    else
        error('Input size mismatch')
    end
    
% ----------------------------------------------------------------------- %
% ELSE
% ----------------------------------------------------------------------- %
else
    error('Cannot multiply these two objects. Only legpoly x scalar (numeric or sdpvar) is supported.')
end


%% END SCRIPT
end