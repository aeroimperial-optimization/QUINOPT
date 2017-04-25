function p = mrdivide(x,y)

%% OVERLOADED: dvarpoly/mrdivide

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    25/04/2017
% Last Modified:    25/04/2017
% ----------------------------------------------------------------------- %

if numel(y)==1
    % Reuse code!
    y = (1/y).*ones(size(x));
    p = times(x,y);
    return
    
else
    error('Cannot divide a dvarpoly object by a matrix. Only division by scalars is allowed.')
    
end
    
    