function display(x)

% OVERLOADED: indvar/display

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    25/04/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% Input parameters
DOMAIN = getDomain(x);
nvars =  size(DOMAIN,1);

% Build composite domain expression
dom = sprintf('[%.6g,%.6g]',DOMAIN(1,1),DOMAIN(1,2));
for i=2:nvars
    dom = sprintf('%sx[%.6g,%.6g]',dom,DOMAIN(i,1),DOMAIN(i,2));
end

% Display
if isValidindvar(x)
    fprintf('Independent variable, domain %s\n',dom);
else
    % Call yalmip display function
    display(sdpvar(x));
end

end