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
expr = u(x,1)^2 - nu*u(x)^2;

% Then, we set the boundary and symmetry conditions on u(x). The periodic 
% boundary conditions is enforced as u(-1)-u(1)=0, while the symmetry condition
% can be enforced using the command "assume()":
bc = [ u(-1)-u(1) ];
assume(u,'odd')

%% Solve in a loop for different SDP approximation degrees
% 
% We compute upper and lower bounds on the optimal nu using QUINOPT's inner and
% outer approximation methods, using different Legendre expansion degrees 
% (OPTIONS.N) in a loop. To speed up the loop and run YALMIP in silent mode, we
% set some YALMIP options:
OPTIONS.YALMIP = sdpsettings('verbose',0,'cachesolvers',1);

% Then we display a header to tabulate the results
fprintf('\nExample 5:\n')
fprintf('================================\n')
fprintf('  N  |  UB/pi^2   |  LB/pi^2   |\n')
fprintf('================================\n')

% Finally we loop:
for N = 1:10
    
    % Compute upper bound (using OPTIONS.method = 'outer' since we maximize nu)
    OPTIONS.method = 'outer';
    OPTIONS.N = N;
    quinopt(expr,bc,-nu,OPTIONS);
    UB = value(nu)/pi^2;
    
    % When N is too small, the outer approximation may be unbounded, implying an
    % infinite upper bound!
    if  isnan(UB); UB = +inf; end

    % Compute lower bound (using OPTIONS.method = 'inner' since we maximize nu).
    % First we reset the value of nu stored by YALMIP to be able to detect
    % unboundedness/infeasibility.
    assign(nu,NaN);                         
    OPTIONS.method = 'inner';
    sol = quinopt(expr,bc,-nu,OPTIONS);
    LB = value(nu)/pi^2;
    if  isnan(LB); LB = -inf; end
    
    % Finally print the results to screen
    fprintf(' %2i  |  %8.6f  |  %8.6f  |\n',N,UB,LB)
    
end
fprintf('================================\n')

%% END CODE