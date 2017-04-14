function clearModel

% CLEARMODEL.m Clear variables for QUINOPT
%
% CLEARMODEL clears all internal variables used by QUINOPT. Moreover,
%       any variables of class <a href="matlab:help('legpoly')">legpoly</a>, <a href="matlab:help('depvar')">depvar</a> and <a href="matlab:help('indvar')">indvar</a> are removed from 
%       the workspace from which CLEARMODEL is called.
%
% NOTE: Consider running the "yalmip clear" after CLEARMODEL to avoid
%       performance issues due to a build-up of unused YALMIP variables.
% 
% See also YALMIP, QIIMODEL, DVARPOLY, DVARLIST

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    05/05/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% Clear variables from workspace
W = evalin('caller','whos');
for i = 1:size(W,1)
    if any(strcmp(W(i).class,{'legpoly','indvar','depvar'}))
        evalin('caller', ['clear ' W(i).name ';']);
    end
end

% Clear persistend variables
evalin('caller','dvarlist clear'); 
evalin('caller','qiimodel clear');