function s = plus(x,y)

%% OVERLOADED: dvarpoly/plus

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    16/04/2015
% Last Modified:    16/04/2016
% ----------------------------------------------------------------------- %

% Check for size mismatch
[mx,nx] = size(x);
[my,ny] = size(y);
if (numel(x)~=1 && numel(y)~=1) && (mx~=my || nx~=ny)
    error('Size mismatch.')
end

% Check to expand scalar input to matrix size
if numel(x)==1 && numel(y)~=1
    x = repmat(x,my,ny);
    
elseif numel(y)==1 && numel(x)~=1
    y = repmat(y,mx,nx);
    
end

% Operation
if isa(x,'dvarpoly') && isa(y,'dvarpoly')
    
    s = struct([]);
    for i=1:mx
        for j=1:nx
            
            % Find monomials
            [ivars,monom_x,monom_y] = findCommonBase(x(i,j),y(i,j));
            
            % Concatenation
            coeff = [x(i,j).coeff; y(i,j).coeff];
            monom = [monom_x; monom_y];
            
            % Find repeated monomials
            [monom,ia,ic] = unique(monom,'rows','first');
            maxind = max(ic);
            newcoeff = [];
            for k = 1:maxind
                rc = coeff(ic==k);
                newcoeff = [newcoeff; sum(rc)];    % works with sdpvars
            end
            
            % Find monomials with non-zero coefficients
            if isnumeric(newcoeff)
                ind = newcoeff~=0;
            else
                % if sdpvar or legpoly
                ind = any(newcoeff);
            end
            
            % Set monomials
            if any(ind)
                s(i,j).coeff = newcoeff(ind);
                s(i,j).ivars = ivars;
                s(i,j).monom = monom(ind,:);
            else
                s(i,j).coeff = 0;
                s(i,j).ivars = [];
                s(i,j).monom = 0;   
            end
            
        end
    end
    s = dvarpoly(s);
    
elseif isa(x,'dvarpoly')
    % added a constant y to a dvarpoly x
    s = struct([]);
    if isa(y,'indvar'); y=sdpvar(y); end
    for i=1:mx
        for j=1:nx
            
            if isZero(y(i,j))
                s(i,j).coeff = x(i,j).coeff;
                s(i,j).ivars = x(i,j).ivars;
                s(i,j).monom = x(i,j).monom;
                
            elseif isZero(x.monom(1,:))
                % add to the coefficient
                s(i,j).coeff = x(i,j).coeff;
                s(i,j).coeff(1) = s(i,j).coeff(1)+y(i,j);
                s(i,j).ivars = x(i,j).ivars;
                s(i,j).monom = x(i,j).monom;
                
            else
                s(i,j).coeff = [y(i,j); x(i,j).coeff];
                s(i,j).ivars = x(i,j).ivars;
                s(i,j).monom = [zeros(1,length(x(i,j).ivars)); x(i,j).monom];
                
            end
            
        end
    end
    s = dvarpoly(s);
    
elseif isa(y,'dvarpoly')
    % Reuse code
    s = plus(y,x);
    return
    
end