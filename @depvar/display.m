function display(u)

% OVERLOADED: depvar/display

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    25/04/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% Extract model for this variable
mod = qiimodel('query');
varInd = find(mod.DEPVARMODEL.DVARID==u.id);
if isempty(varInd), error('Invalid dependent variable identifier.'); end
MAXDER = mod.DEPVARMODEL.MAXDER(varInd);
fprintf('Dependent variable (%i derivatives defined)\n',MAXDER);