function p = times(x,y)

%% OVERLOADED: dvarpoly/times
% Element-wise multiplication

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
    
    % Loop over entries to operate - quicker than arrayfun?
    p = struct([]);
    for i=1:mx
        for j=1:nx
            
            if isZero(x(i,j).coeff) || isZero(y(i,j).coeff)
                p(i,j).coeff = 0;
                p(i,j).ivars = [];
                p(i,j).monom = 0;
                
            else
                % Find coefficients (outer product if vectors)
                %coeff = x(i,j).coeff * y(i,j).coeff'; WRONG?
                coeff = y(i,j).coeff * x(i,j).coeff';
                
                % Find monomials
                [ivars,monom_x,monom_y] = findCommonBase(x(i,j),y(i,j));
                nmonsx = size(monom_x,1);
                nmonsy = size(monom_y,1);
                indx = repmat(1:nmonsx,nmonsy,1);
                monom = monom_x(indx(:),:) + repmat(monom_y,nmonsx,1);
                
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
                    p(i,j).coeff = newcoeff(ind);
                    p(i,j).ivars = ivars;
                    p(i,j).monom = monom(ind,:);
                else
                    p(i,j).coeff = 0;
                    p(i,j).ivars = [];
                    p(i,j).monom = 0;
                end
                
                
            end
        end
    end
    
    
elseif isa(x,'dvarpoly')
    % Only need to multiply coefficients!
    % Loop over entries to operate - quicker than arrayfun?
    p = struct([]);
    if isa(y,'indvar'); y=sdpvar(y); end
    for i=1:mx
        for j=1:nx
            
            if isZero(y(i,j))
                % Set to zero!
                p(i,j).coeff = 0;
                p(i,j).ivars = [];
                p(i,j).monom = 0;
            else
                % Multiply coefficients only
                p(i,j).coeff = x(i,j).coeff * y(i,j);
                p(i,j).ivars = x(i,j).ivars;
                p(i,j).monom = x(i,j).monom;
            end
            
        end
    end
    
    
elseif isa(y,'dvarpoly')
    % Reuse code
    p = times(y,x);
    return
    
    
end

% Assign class
p = dvarpoly(p);

% END FUNTION
end