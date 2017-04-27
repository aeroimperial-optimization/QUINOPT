%% example04.m
%
% Find the minimum constant gamma for which the inequality
%
% /1
% |
% |  gamma*u_xx^2 - u_x^2 dx  >= 0
% |
% /-1
%
% holds for all functions u(x) subject to the boundary conditions u(-1)=0,
% u(1)=0, u_x(-1)-u_x(1)=0.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    24/02/2016
% Last Modified:    27/04/2017
% ----------------------------------------------------------------------- %

%% CODE

clear;             % clean workspace 
yalmip clear;      % clear YALMIP's internal variables
quinopt clear;     % clear QUINOPT's internal variables

%% Initialization
x = indvar(-1,1); % Initialize dependent variable with domain [-1,1]
u = depvar(x);    % Dependent on x
parameters gamma

%% Problem setup
expr = gamma*u(x,2)^2 - u(x,1)^2;
BC = [u(-1); u(1); u(-1,1)-u(1,1)];

%% Solve optimization to minimize nu
% We solve the problem for a number of values of theLegendre series 
% truncation parameter N, starting with N=2 (the minimum). 

% To run in silent mode, we set YALMIP's option 'verbose' to 0. Also, we
% speed up the iteration by settin YALMIP's 'cachesolvers' option to 1.
opts.YALMIP = sdpsettings('verbose',0,'cachesolvers',1);

% Let us display a nice header...
fprintf('\n==============================================================\n')
fprintf('|   N   |   gamma_opt   |   0.25*pi^2-gamma_opt   |   time   |\n')
fprintf('==============================================================\n')
for N = 2:7
    
    % Call and time QUINOPT
    time = tic;
    opts.N = N;
    quinopt(expr,BC,gamma,opts);
    time = toc(time);

    % Extract the solution and compare it to the analytical answer
    gamma_opt = value(gamma);
    error = 1/pi^2 - gamma_opt;
    
    % Print information
    fprintf('|   %i   |   %8.6f    |       %11.4e       |  %6.4f  |\n',N,gamma_opt,error,time)
end
fprintf('==============================================================\n\n')

%% END SCRIPT