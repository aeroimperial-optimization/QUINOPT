%% AdditionalEx1.m
%
% Find the minimum constant nu for which the inequality
%
% /1
% |
% |  u_xx^2 - nu*u_x^2 dx  >= 0
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
% Last Modified:    12/05/2016
% ----------------------------------------------------------------------- %

%% CODE

clear;          % clean workspace 
clearModel;     % clear QUINOPT's internal variables

%% Initialization
x = indvar(-1,1); % Initialize dependent variable with domain [-1,1]
u = depvar(x);    % Dependent on x
parameters nu

%% Problem setup
expr = nu*u(x,2)^2 - u(x,1)^2;
BC = [u(-1); u(1); u(-1,1)-u(1,1)];

%% Solve optimization to minimize nu
quinopt(expr,BC,nu);
nuval = value(nu);
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('AdditionalEx1: DEFAULT SOLUTION INFO');
disp(['Numerical Solution       :        nu = ',num2str(nuval)])
disp(['Error with analytical ans: nu-1/pi^2 = ',num2str(nuval-1/pi^2)])
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++'); disp(' ')

%% Refine solution
quinopt(expr,BC,nu,[],[],5);
nuval = value(nu);
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('AdditionalEx1: REFINED SOLUTION INFO');
disp(['Numerical Solution       :        nu = ',num2str(nuval)])
disp(['Error with analytical ans: nu-1/pi^2 = ',num2str(nuval-1/pi^2)])
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++'); disp(' ')

