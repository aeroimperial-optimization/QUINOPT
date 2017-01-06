function [C,IA,IC] = getVariables(A,flag)

%% GETVARIABLES.m Returns variables in dvarpoly
%
% [C,IA,IC] = GETVARIABLES(A,flag) returns indices to the variables that
%       appear in the dvarpoly A. By default, the variable indices are
%       sorted. To preserve the original order of variables in A, specify
%       'stable' as second input.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    16/04/2015
% Last Modified:    06/01/2017
% ----------------------------------------------------------------------- %

% Compatible version with old MATLAB versions
if nargin < 2
    % flag = 'sorted';
    %[C,IA,IC] = unique([A.ivars],'first');
    [C,IA,IC] = unique(vertcat(A.ivars).','first');
    
elseif strcmpi(flag,'stable')
%     x = [A.ivars];
%     [C,IA,IC] = unique(x,'first');
    x = vertcat(A.ivars);
    [C,IA,IC] = unique(x,'first');
    C = x(sort(IA)).';
    
elseif strcmpi(flag,'rows')
     error('Flag "rows" not allowed.')
     
end