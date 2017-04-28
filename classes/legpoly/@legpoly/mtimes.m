function s = mtimes(X,Y)

%% mtimes.m
%
% c = mtimes(X,Y)

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

%% CODE

if numel(X)==1 || numel(Y)==1
    % If either factor is a scalar
    s = X.*Y;
else
    [rowsX,colsX] = size(X);
    [rowsY,colsY] = size(Y);
    if colsX~=rowsY
        error('Dimension of factors mismatch.')
    end
    s = [];
    X = X.';
    for j=1:colsY
        for i=1:rowsX
           s = [s; sum( X(:,i).*Y(:,j) )]; 
        end
    end
    s = reshape(s,rowsX,colsY);
    
%% END SCRIPT
end