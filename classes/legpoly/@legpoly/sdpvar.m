function Y = sdpvar(X)

% OVERLOADED: legpoly/sdpvar
%
% Convert legendre polynomial object into an sdpvar polynomial

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

[m,n] = size(X);
Y = [];

for i = 1:numel(X)
    
    z = X(i).ivar;                     % get independent variable in [-1,1]
    DOMAIN = X(i).domn;                % get domain
    
    if any(DOMAIN)
        % not a constant
        c = monBasisCoef(X(i).coef);       % find coefficients for polynomial over [-1,1]
        deg = length(c)-1;
        mons = monolist(z,deg,0);
        p = mons'*c;
        
        % Rescaled variables
        x = (2*z-DOMAIN(2)-DOMAIN(1))/(DOMAIN(2)-DOMAIN(1));
        p = replace(p,z,x);
        
        % add polynomial to list
        Y = [Y; p];  
        
    else
        % a constant
        Y = [Y; X(i).coef];
    end
    
end
Y = reshape(Y,m,n);