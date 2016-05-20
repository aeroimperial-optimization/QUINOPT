function [xInList,ListIndex] = belongsTo(x,list)

% Check if vector of sdpvar variables "list" contains any of the sdpvar
% variables listed in the vector "x". Returns:
%
% - xInList: a logic vector with 1 corresponding to elements of x which are in
%            list, 0 corresponding to those that are not listed in list
% - ListIndex: a vector with the first element in list which is equal to the
%              corresponding element in x. For example, ListIndex = [0 4 2]
%              means that x(1) is not in the list, x(2)==list(4) and
%              x(3)==list(2).

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    16/03/2015
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% Vectorise inputs
x = x(:);
list = list(:);

% Find basis and variables in x
xVars = depends(x);
xBase = getbase(x);

% Find basis and variables in list
lVars = depends(list);
lBase = getbase(list);

% Find list of all variables
allVars = sort(unique([xVars,lVars]));

% Inflate xBase to include all variables
rows = size(xBase,1);
cols = length(allVars)+1;
[tmp,p] = ismember(xVars,allVars);
[n,m] = meshgrid([1; p'+1],(1:rows)');
xBase = sparse(m(:),n(:),xBase,rows,cols);

% Inflate lBase to include all variables
rows = size(lBase,1);
[tmp,p] = ismember(lVars,allVars);
[n,m] = meshgrid([1; p'+1],(1:rows)');
lBase = sparse(m(:),n(:),lBase,rows,cols);

% Find if any rows of xBase are rows of lBase
if nargout==1
    [xInList] = ismember(xBase,lBase,'rows');
elseif nargout==2
    [xInList,ListIndex] = ismember(xBase,lBase,'rows');
else
    error('Too many outputs.')
end

