%% example05.m
%
% Find the minimum constant nu for which the Poincare inequality
%
% /1
% |
% |  u_x^2 - nu*u^2 dx  >= 0
% |
% /-1
%
% holds for all functions u(x) subject that are odd and satisfy the periodicity
% condition u(-1)=u(1). The analytical answer is nu = pi^2.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    10/05/2017
% Last Modified:    10/05/2017
% ----------------------------------------------------------------------- %

%% CODE

%% Create the variables
%
% As usual, we start by cleaning the workspace and the internal variables in
% QUINOPT and YALMIP:
clear
yalmip clear
quinopt clear

% Then we create the integration variable x, the dependent variable u(x), and
% the optimization variable nu:
x = indvar(-1,1);
u = depvar(x);
parameters nu

%% Set up the inequality
%
% The first step to define the inequality constraint is to construct the
% integrand of the inequality:
EXPR = u(x,1)^2 - nu*u(x)^2;

% Then, we set the boundary and symmetry conditions on u(x). The periodic
% boundary conditions is enforced as u(-1)-u(1)=0, while the symmetry condition
% can be enforced using the command "assume()":
BC = [ u(-1)-u(1) ];
assume(u,'odd')

%% Solve the problem
%
% To solve the problem and maximize nu, we use the command "quinopt()" with three
% arguments: EXPR, BC and the objective function. Since QUINOPT minimizes the
% specified objective function, we maximize nu by minimizing -nu.
quinopt(EXPR,BC,-nu);
value(nu)/pi^2

% We can also refine the approximation of the integral inequalities by increasing
% the number of Legendre coefficients used by QUINOPT to expand them. We do this
% by setting the option "N", as in example04.m:
options.N = 5;
quinopt(EXPR,BC,-nu,options);
value(nu)/pi^2

%% END CODE