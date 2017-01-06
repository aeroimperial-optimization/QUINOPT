%% AdditionalEx4.m
%
% Compute the optimal background field for 3D stress-driven shear flow at
% Gr=1000. We consider perturbations at the first 3 wave numbers. This 
% example illustrates how to solve a problem with multiple integral 
% inequality constraints.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/04/2016
% Last Modified:    12/05/2016
% ----------------------------------------------------------------------- %

%% CODE

clear; clearModel;

%% Parameters
Gamma = 2;   % period in x direction, x in [0,xPeriod]
P = 20;      % degree of a polynomial background field.
G = 1000;    % value of forcing parameter


%% Setup problem variables
x = indvar(-1,1);                   % the independent variable
[u,v] = depvar(x);                  % 2 dependent variables

%% Setup Phi_x as a sum of Legendre polynomials
% Phi_x : Legendre polynomial representing the derivative of the background field
% PhiHat: the Legendre coefficients of Phi_x 
%
% NOTE: the optimization variables used in this problem are the Legendre
% coefficients PhiHat, which implicitly define the background field Phi_x
[Phi_x,PhiHat] = legpoly(x,P-1);

%% Setup inequality
expr = [];
for n=1:3
    k = 2*pi*n/Gamma;   % the wave number
    expr = [expr; 16/k^2*v(x,2)^2 + 8*v(x,1)^2 + 4*u(x,1)^2 + k^2*( u(x)^2+v(x)^2 ) ...
            + 4*Phi_x/k*u(x)*v(x)];
end
BC = [u(-1); u(1,1)];                       % BC on u
BC = [BC; v(-1); v(1); v(-1,1); v(1,2)];    % BC on v  

%% Setup objective and additional constraints
parameters eta
obj = 2*eta - 4*PhiHat(1);                                % the objective
S = [diag( (2*(0:P-1)+1)/2 ), PhiHat; PhiHat', G*eta];    % matrix for additional constraints
cnstr = [S>=0; sum(PhiHat)==G/2];                         % the additional constraints

%% Solve optimization with integral inequality and additional constraint
time = tic;
quinopt(expr,BC,obj,cnstr);
time = toc(time);

%% Compute and plot background field
z = -1:0.001:1;
Phi_x = replace(Phi_x,PhiHat,value(PhiHat));        % assign value to coefficients
Phi = legpolyint(Phi_x,x,0);                        % integrate Phi_x with BC Phi(-1)=0
PhiVal = legpolyval(Phi,z);                         % evaluate Phi at z

figure(gcf); clf; plot(PhiVal/PhiVal(end),z,'-','LineWidth',1.5); 
axis equal; xlim([0,1]); ylim([-1,1]); 
xlabel('\phi(z)/\phi(1)'); ylabel('z'); title('Example 5');

%% Display results
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('QUINOPT example 4: SOLUTION INFO');
disp(['Solution time:  ',num2str(time),' seconds'])
disp(['Optimal objective value = ',num2str(value(obj))])
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++'); disp(' ')