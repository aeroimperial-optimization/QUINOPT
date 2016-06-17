function [allvars,monom_x,monom_y] = findCommonBase(x,y)

%% FINDCOMMONBASE.m Common basis for two dvarpolys

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    16/04/2015
% Last Modified:    16/04/2016
% ----------------------------------------------------------------------- %

if numel(x)~=1 && numel(y)~=1
    error('Function "findcommonbase" only works for scalar inputs.')
end

% Get independent variables &  union
ivars_x = x.ivars;
ivars_y = y.ivars;
allvars = union(ivars_x,ivars_y);       % already sorted
nvars = length(allvars);

% Inflate monomials
num_monom_x = size(x.monom,1);
num_monom_y = size(y.monom,1);
[tmp,indx] = ismember(ivars_x,allvars);
[tmp,indy] = ismember(ivars_y,allvars);
monom_x = zeros(num_monom_x,nvars);
monom_x(:,indx) = x.monom;
monom_y = zeros(num_monom_y,nvars);
monom_y(:,indy) = y.monom;