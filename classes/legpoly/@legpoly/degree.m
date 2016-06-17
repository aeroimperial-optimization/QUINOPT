function DEG = degree(p,x,flag)

%% DEGREE.m
%
% DEG = degree(p,x)

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

%% CODE
if nargin < 2
    x = 0;
    flag = 'max';
elseif nargin < 3
    flag = 'max';
end

DEG = cellfun(@(z)length(z)-1,{p.coef}');
if strcmpi(flag,'all')
    [m,n] = size(p);
    DEG = reshape(DEG,m,n);
    return
else
    DEG = max(DEG(:));
end


end
%% END SCRIPT
