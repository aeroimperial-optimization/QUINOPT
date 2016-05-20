function [coeff,monom,ivars] = coefficients(x)

%% OVERLOADED: dvarpoly/coefficients
% Return all coefficients, monmials and variables of the dvarpoly x

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    16/04/2015
% Last Modified:    16/04/2016
% ----------------------------------------------------------------------- %

if numel(x) > 1
    [m,n] = size(x);
    [coeff{1:m,1:n}] = deal(x.coeff);
    [monom{1:m,1:n}] = deal(x.monom);
    [ivars{1:m,1:n}] = deal(x.ivars);
    
else
    coeff = x.coeff;
    monom = x.monom;
    ivars = x.ivars;
end