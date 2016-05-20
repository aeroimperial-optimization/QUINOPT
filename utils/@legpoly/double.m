function Ld = double(L)

%% Convert to double - give NaN by default

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

if degree(L)==0
    Ld = coefficients(L);
else
    [m,n]=size(L);
    Ld = NaN*ones(m,n);
end