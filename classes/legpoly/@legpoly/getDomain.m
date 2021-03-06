function DOMAIN = getDomain(P)

%% GETDOMAIN.m Get domain of legpoly
%
% DOMAIN = GETDOMAIN(P) returns the domain of the independent variable of a 
%   polynomial in the Legendre basis (class <a href="matlab:help('legpoly')">legpoly</a>). The domain [0,0] is 
%   reserved to constants (polynomials of degree 0).
%
% See also INDVAR, LEGPOLY

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/04/2016
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

ivar = getivar(P); % Find the independent variables of each entry
[i,j]=find(any(ivar),1,'first');   % indices of polynomial of degree >=1
if ~isempty(i)
    DOMAIN = P(i,j).domn;
else
    DOMAIN = [0,0];
end