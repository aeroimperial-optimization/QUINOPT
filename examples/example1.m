%% Example1.m
%
% Find the minimum constant nu for which the Poincare inequality
%
% /1
% |
% |  u_x^2 - nu*u^2 dx  >= 0
% |
% /-1
%
% holds for all functions u(x) subject to the Dirichlet boundary conditions 
% u(-1)=0, u(1)=0.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    24/02/2016
% Last Modified:    02/05/2016
% ----------------------------------------------------------------------- %


%% Initialization
%
% First, we remove any existing variables to avoid unexpected dependencies
clear;          % clear workspace
clearModel;     % clear QUINOPT's internal variables

% Then we initialize the independent variable of integration, the dependent
% variable u and the optimization variable gamma with the commands
x = indvar(-1,1);  % define independent variable with domain [-1,1]
u = depvar(x);     % define u, which depends on x
parameters gamma   % define the optimsation variable gamma

%% Construct the problem
% The integrand of the integral inequality is constructed with the commands
expr = u(x,1)^2 - gamma*u(x)^2;
%
% where u(x,1) specifies the first derivative of u evaluated at the
% variable point x. 

% The boundary conditions u(-1)=0 and u(1)=0 are specified by constructing 
% the vector
BC = [u(-1); u(1)];
%
% where the syntax u(0) specifies the value of $u$ at the boundary point 
% x = -1.

% Finally we specify the objective to minimize as
obj = -gamma;

%% Solve the problem 
% We solve the problem for a number of values of theLegendre series 
% truncation parameter N. 

%To run in silent mode, we set YALMIP's option 'verbose' to 0. Also, we
%speed up the iteration by settin YALMIP's 'cachesolvers' option to 1.
opts.YALMIP = sdpsettings('verbose',0,'cachesolvers',1);

% Display a nice header
fprintf('|============================================================|\n')
fprintf('|   N   |   gamma_opt   |   0.25*pi^2-gamma_opt   |   time   |\n')
fprintf('|============================================================|\n')
for N = 1:7
    
    % Call and time QUINOPT
    time = tic;
    quinopt(expr,BC,obj,[],[],N,opts);
    time = toc(time);

    % Extract the solution and compare it to the analytical answer
    gamma_opt = value(gamma);
    error = 0.25*pi^2 - gamma_opt;
    
    % Print information
    fprintf('|   %i   |   %8.6f    |       %11.4e       |  %6.4f  |\n',N,gamma_opt,error,time)
end
fprintf('|============================================================|\n')

%% END SCRIPT
