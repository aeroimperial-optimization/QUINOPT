function param = findSosParameters(p,x)

%% FINDSOSPARAMETERS.m Find parameters in polynomial.
%
% param = FINDSOSPARAMETERS(p,x) returns a vector containing the YALMIP
%       variables that appear in the polynomial p in addition to the
%       polynomial independent variables x (x can be a vector of multiple
%       independent variables)

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    11/04/2016
% Last Modified:    11/04/2016
% ----------------------------------------------------------------------- %

% Find variables and base
c = coefficients(p,x);                          % coefficients of p
cvars = depends(c);                             % find ID of YALMIP variables
nvars = length(cvars);                          % number of variables
B = spdiags(ones(nvars,1),1,nvars,nvars+1);     % base matrix

% Create YALMIP variable using YALMIP internal calling syntax
param = sdpvar(nvars,1,[],cvars,B);             