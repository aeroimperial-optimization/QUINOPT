%% example04.m
%
% Find the minimum constant C for which the inequality
%
% /2*pi
% |
% |  C*u_xx^2 - u_x^2 dx  >= 0
% |
% /0
%
% holds for all functions u(x) subject to the boundary conditions u(0)=0,
% u(2*pi)=0, u_x(0)-u_x(2*pi)=0.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    24/02/2016
% Last Modified:    27/04/2017
% ----------------------------------------------------------------------- %

%% CODE

%% Initialization
%
% First, we remove any existing variables to avoid unexpected dependencies
clear;              % clear workspace
yalmip clear;       % clear YALMIP's internal variables
quinopt clear;      % clear QUINOPT's internal variables

% Then we initialize the variables:
x = indvar(0,2*pi);        % Initialize integration variable with domain [0,2*pi]
u = depvar(x);             % The dependent variable
parameters C               % The optimization parameter

%% Setting up the integral inequality
%
% The integrand of the integral inequality and the boundary conditions on u are
% set up with the commands
expr = C*u(x,2)^2 - u(x,1)^2;
bc = [u(0); u(2*pi); u(0,1)-u(2*pi,1)];


%% Optimize
% We seek upper and lower bounds on the lowest possible C. QUINOPT's default 
% behaviour is to compute upper bounds on the optimal objective value by 
% formulating and inner approximation of the feasible set of the integral 
% inequality constraints. A lower bound is found by overriding the default
% method with an "options" structure defined as
options.method = 'outer';
% This makes QUINOPT formulate an outer approximation of the feasible set of the
% integral inequality constraints. 

% We also tell QUINOPT to run quietly by setting YALMIP's "verbose" option to 0,
% and cache YALMIP's available solvers for speed
options.YALMIP = sdpsettings('verbose',0,'cachesolvers',1);

% We can then solve the problem:
[sol,cnstr,data] = quinopt(expr,bc,C,options);
LB = value(C);                      % extract the lower bound on the optimal C

% To compute an upper bound, we need to reset QUINOPT's default behaviour:
options.method = 'inner';
quinopt(expr,bc,C,options);
UB = value(C);                      % extract the upper bound on the optimal C

%% Display
fprintf('\nExample 4:\n')
fprintf('=================================================================\n')
fprintf('|        N        |       LB      |      UB      |  GAP (UB-LB) |\n')
fprintf('=================================================================\n')
fprintf('|   2 (default)   |   %8.6f    |   %8.6f   |   %8.6f   |\n',LB,UB,UB-LB)


%% Improve the bounds
% To improve the bounds, we can refine the approximation of the integral
% inequalities by increasing the number of Legendre coefficients used by QUINOPT
% to expand them. We do this by setting the option "N":
for N = 3:9
    
    % Set number of Legendre coefficients
    options.N = N;
    
    % Compute lower bound
    options.method = 'outer';
    quinopt(expr,bc,C,options);
    LB = value(C);
    
    % Compute upper bound
    options.method = 'inner';
    quinopt(expr,bc,C,options);
    UB = value(C);
    
    % Print information
    fprintf('|      %4i       |   %8.6f    |   %8.6f   |   %4.2e   |\n',N,LB,UB,UB-LB)
end
fprintf('=================================================================\n')

%% END CODE