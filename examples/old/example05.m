%% example3.m
%
% Compute the optimal background field for 2D stress-driven shear flow at
% Gr=1000 considering only one wavenumber k. This example illustrates how 
% to specify constraints on the optimization variables in addition to the
% integral inequality constraint.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/04/2016
% Last Modified:    17/04/2017
% ----------------------------------------------------------------------- %

%% CODE

clear;             % clean workspace 
yalmip clear;      % clear YALMIP's internal variables
quinopt clear;     % clear QUINOPT's internal variables

%% Parameters
Gamma = 2;   % period in x direction, x in [0,xPeriod]
P = 20;      % degree of polynomial Phi
G = 1000;    
k = 2*pi/Gamma;  % the first wavenumber


%% Setup problem variables
x = indvar(-1,1);               % the independent variable
[u,v] = depvar(x);              % the dependent variables

%% Setup Phi_x as a sum of Legendre polynomials
% Phi_x: the actual polynomial, represented in the Legendre basis
% PhiHat: the Legendre coefficients of Phi
%
% NOTE: the optimization variables used in this problem are the Legendre
% coefficients PhiHat, which implicitly define the background field Phi_x
[Phi_x,PhiHat] = legpoly(x,P-1);

%% Setup integral inequality constraint
expr = 16/k^2*( u(x,2)^2+v(x,2)^2 ) + 8*( u(x,1)^2+v(x,1)^2 ) + ...
       k^2*( u(x)^2+v(x)^2 )  - 8*Phi_x/k*( u(x)*v(x,1) - u(x,1)*v(x) );
BC = [u(-1); u(1); u(-1,1); u(1,2)];        % BC on u
BC = [BC; v(-1); v(1); v(-1,1); v(1,2)];    % BC on v  

%% Setup objective and additional constraints
parameters eta                                            % the additional problem variable
obj = 2*eta - 4*PhiHat(1);                                % the objective function
S = [diag( (2*(0:P-1)+1)/2 ), PhiHat; PhiHat', G*eta];    % matrix for additional constraint
cnstr = [sum(PhiHat)==G/2, S>=0];                         % vector of additional constraint

%% Solve optimization with integral inequality and additional constraint
time = tic;
quinopt(expr,BC,obj,[],cnstr);
time = toc(time);

%% Compute and plot integral of Phi_x 
z = -1:0.001:1;                                 % coordinate values for plotting
Phi_x = replace(Phi_x,PhiHat,value(PhiHat));    % assign value to coefficients
Phi = legpolyint(Phi_x,x,0);                    % integrate Phi with BC Phi(-1)=0
PhiVal = legpolyval(Phi,z);                     % evaluate Phi at z

figure(); clf; plot(PhiVal/PhiVal(end),z,'-','LineWidth',1.5); 
axis equal; xlim([0,1]); ylim([-1,1]); 
xlabel('\phi(x)/\phi(1)'); ylabel('x'); title('Example 4');

%% Display results
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('QUINOPT example 3: SOLUTION INFO');
disp(['Solution time:  ',num2str(time),' seconds'])
disp(['Optimal objective value = ',num2str(value(obj))])
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++'); disp(' ')

  