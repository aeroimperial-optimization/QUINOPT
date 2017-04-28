%% example03.m
%
% Find the minimum constant gamma for which the Poincare inequality
%
% /1
% |
% |  u_x^2 - gamma*u^2 dx  >= 0
% |
% /0
%
% holds for all functions u(x) subject to the Dirichlet boundary conditions 
% u(0)=0, u(1)=0. The analytical answer is gamma = pi^2.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    24/02/2016
% Last Modified:    27/04/2017
% ----------------------------------------------------------------------- %


%% Initialization
%
% First, we remove any existing variables to avoid unexpected dependencies
clear;              % clear workspace
yalmip clear;       % clear YALMIP's internal variables
quinopt clear;      % clear QUINOPT's internal variables

% Then we initialize the independent variable of integration, the dependent
% variable u and the optimization variable gamma with the commands
x = indvar(0,1);       % define independent variable with domain [-1,1]
u = depvar(x);          % define u, which depends on x
parameters gamma        % define the optimsation variable gamma

%% Construct the problem
% The integrand of the integral inequality is constructed with the commands
expr = u(x,1)^2 - gamma*u(x)^2;

% where u(x,1) specifies the first derivative of u evaluated at the
% variable point x. 

% The boundary conditions u(-1)=0 and u(1)=0 are specified by constructing 
% the vector
BC = [u(0); u(1)];

% where the syntax u(0) specifies the value of u at the boundary point x = 0.

% Finally we specify the objective. QUINOPT minimizes the objective function it
% receives, so to maximize gamma we minimize:
obj = -gamma;

%% Solve the problem 
% To minimize the objective function, we use the command quinopt() with three
% inputs: the integrand of the inequality constraint, the boundary conditions,
% and the objective function.
quinopt(expr,BC,obj);
gamma_opt = value(gamma);
error = gamma_opt - pi^2;
    
%% Display results
fprintf('\n======================================\n')
fprintf('Analytical optimum  : gamma = %f\n',pi^2)
fprintf('Optimum from QUINOPT: gamma = %f\n',gamma_opt)
fprintf('Percentage error    : %6.4f %%\n',error/pi^2*100)
fprintf('======================================\n\n')

%% END CODE

%% END SCRIPT