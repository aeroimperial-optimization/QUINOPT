function xm = uminus(x)

%% minus.m
%
% xm = uminus(x)

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
%   Description:    Compute A-B
% ----------------------------------------------------------------------- %

%% CODE

[m,n] = size(x);

[xm(1:m,1:n).coef] = deal(x.coef);
newcoeff = arrayfun(@(x)-x.coef,xm,'uniformoutput',0);
[xm(1:m,1:n).coef] = deal(newcoeff{:});

[xm(1:m,1:n).ivar] = deal(x.ivar);
[xm(1:m,1:n).domn] = deal(x.domn);


xm = legpoly(xm);

end
%% END SCRIPT
