function display(p)

%% DISPLAY (overloaded)
% Display legendre polynomial object

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

[m,n] = size(p); 

% Get number of sdpvar variables in the Legendre coefficients
vars = depends(p(1).coef);
for i = 2:numel(p)
    vars = unique([vars, depends(p(i).coef)]);
end

% Get degree and number of variables
deg = degree(p);
nvars = length(vars);
dom = getDomain(p);

% Display
fprintf('Polynomial in Legendre basis %ix%i (degree %i, %i variables, domain [%.6g %.6g])\n',...
    m,n,deg,nvars,dom(1),dom(2));