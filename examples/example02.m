%% example02.m
%
% Find a feasible value gamma such that
%
% /1
% |
% |  u_x^2 + gamma*u*u_x + u^2 dx  >= 0
% |
% /-1
%
% holds for all functions u(x) that satisfy the boundary conditions u'(0)=0 and
% u(1)=0.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    27/04/2017
% Last Modified:    27/04/2017
% ----------------------------------------------------------------------- %


%% Initialization
%
% First, we remove any existing variables to avoid unexpected dependencies
clear;              % clear workspace
yalmip clear;       % clear YALMIP's internal variables
quinopt clear;      % clear QUINOPT's internal variables

% Then we initialize the independent variable of integration, the dependent
% variable u and the optimization variable gamma:
x = indvar(0,1);       % define independent variable with domain [-1,1]
u = depvar(x);          % define u, which depends on x
parameters gamma        % define the optimzation variable gamma

%% Setting up the inequality
% The integrand of the integral inequality is constructed with the commands
EXPR = u(x,1)^2 +gamma*u(x,1)*u(x) + u(x)^2;

% where u(x,1) specifies the first derivative of u evaluated at the
% variable point x.

% The boundary conditions are specified in a vector BC, such that each boundary
% condition is given by BC(i)=0:
BC = [u(1), u(0,1)];

% where u(0,1) specifies the first derivative of u evaluated at the
% boundary point x=0.

%% Computing a feasible point with QUINOPT
% To compute a feasible value for gamma we use the command quinopt(), and
% YALMIP's command value(). This time, we pass two inputs to quinopt(), the
% integrand of the inequality and the boundary conditions:
quinopt(EXPR,BC);                            % Solve the problem
gamma_feas = value(gamma);                   % Extract the value of gamma

%% Display
fprintf('\nExample 2:\n')
fprintf('=================================\n')
fprintf('Feasible value: gamma = %g\n',gamma_feas)
fprintf('=================================\n')

%% END CODE