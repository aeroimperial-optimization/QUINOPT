%% example08.m
%
% In this example, we compute an upper bound on the energy dissipation for a 2D 
% flow driven by a non-dinemsional surface stress of magnitude G=1000 using the 
% background field method with only one wavenumber k. 
%
% This example illustrates how to specify constraints on the optimization 
% variables in addition to integral inequality constraints.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/04/2016
% Last Modified:    08/05/2017
% ----------------------------------------------------------------------- %

%% CODE
%% Initialization
% First, we clean up
clear;          
yalmip clear;   
quinopt clear;

% Then, we set some problem variables:the horizontal period lambda, the degree
% of the poynomial background field, the Grashoff number G, and the value of the
% first wavenumber.
lambda = 2;  
deg_phi = 20;       
G = 1000;    
k = 2*pi/lambda;


%% Set up the problem variables
% The independent and dependent variables are set up with the usual commands
y = indvar(0,1);               
[u,v] = depvar(y);              

% Then, we construct the polynomial background field phi. Since QUINOPT
% represents polynomials in the Legendre basis internally, we define phi in the
% Legendre basis directly using the command "legpoly()":
phi = legpoly(y,deg_phi);

% To set up the integral inequality constraint, moreover, we need the first
% derivative of phi. This is easily found using the function "jacobian()":
D1phi = jacobian(phi,y);

% Finally, we need to specify the boundary conditions on phi, that are,
% \phi(0)=0 and \phi'(1)=G. The boundary values of phi and its derivatives can
% be computed using the function "legpolyval()":
BC_phi = [legpolyval(phi,0)==0; ...                     % \phi(0)=0
          legpolyval(D1phi,1)==G];                      % \phi'(1)=G

%% Set up the integral inequality constraint
% The integral inequality constraint is set up by specifying the integrand and
% the boundary conditions on the dependent variables u and v:
EXPR = ( u(y,2)^2+v(y,2)^2 )/k^2 + 2*( u(y,1)^2+v(y,1)^2 ) + ...
       k^2*( u(y)^2+v(y)^2 )  - 2*D1phi/k*( u(y)*v(y,1) - u(y,1)*v(y) );
BC = [u(0); u(1); u(0,1); u(1,2)];        % BC on u
BC = [BC; v(0); v(1); v(0,1); v(1,2)];    % BC on v  

%% Set up the objective function
% To set up the objective function, we can use the command "legpolyval()" to
% compute the boundary value \phi(1), and the command "int()" to compute the
% integral of the square of \phi'(z):
OBJ = 2*legpolyval(phi,1) - int(D1phi^2,y,0,1)/G;

%% Solve optimization with integral inequality and additional constraint
time = tic;
quinopt(EXPR,BC,-OBJ,[],BC_phi);
UB = G/(value(OBJ))^2;
time = toc(time);

%% Plot the optimal background field
% Plot the background field and its derivarive using the command "plot()":
subplot(2,1,1)
plot(0:0.01:1,phi,'-','LineWidth',1.5); 
xlabel('$y$','interpreter','latex');
ylabel('$\phi(y)$','interpreter','latex'); 
subplot(2,1,2)
plot(0:0.01:1,D1phi,'-','LineWidth',1.5); 
xlabel('$y$','interpreter','latex');
ylabel('$\phi''(y)$','interpreter','latex');


%% Display results
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('QUINOPT example 3: SOLUTION INFO');
disp(['Solution time:  ',num2str(time),' seconds'])
disp(['Upper bound on dissipation coefficient = ',num2str(UB)])
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++'); disp(' ')

%% END CODE 