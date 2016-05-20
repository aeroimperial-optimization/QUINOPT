function TF = isValidindvar(x)

%% ISVALID.m Test if input is valid independent variable

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    25/04/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

TF = false;
if numel(x)==1
    b = getbase(x);
    d = degree(x);
    if sum(b)==1 && d==1
        TF = true;
    end
end


end