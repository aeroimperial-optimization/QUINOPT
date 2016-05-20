function s = mrdivide(X,Y)

%% mrdivide.m
%
% s = mrdivide(X,Y)  divides a legpoly by a scalar.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

%% CODE

if numel(Y)==1
    % If divisor is a scalar
    s = X./Y;
else
    error('Right matrix division not supported for legpoly objects.')
end
    

%% END SCRIPT
end