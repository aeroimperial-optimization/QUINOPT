function p = mtimes(x,y)

%% OVERLOADED: dvarpoly/mtimes

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    16/04/2015
% Last Modified:    16/04/2016
% ----------------------------------------------------------------------- %

if numel(x)==1 || numel(y)==1
    % Reuse code!
    p = times(x,y);
    return
    
elseif size(x,2)~=size(y,1)
    error('Dimensions of factors mismatch.')
    
else
    % Inefficient loop but simple
    [rowsX,colsX] = size(x);
    [rowsY,colsY] = size(y);
    if colsX~=rowsY
        error('Dimension of factors mismatch.')
    end
    p = [];
    for j=1:colsY
        for i=1:rowsX
           p = [p; sum( x(i,:)'.*y(:,j) )]; 
        end
    end
    p = reshape(p,rowsX,colsY);
    
end
    
    