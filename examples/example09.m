%% example09.m
%
% Compute bounds on energy dissipation for 3D plane Couette flow using the
% indefinite storage functional method (this is equivalent to the usual
% background method formulation of the problem, see e.g. Plasting & Kerswell, J.
% Fluid Mech. 477, 363–379, 2003).
%
% A full description of the problem can be found online at
% http://quinopt.readthedocs.io/05_advanced_examples/planeCouetteBF.html

%% Initialization
%
% First, we remove any existing variables to avoid unexpected dependencies, and
% clear YALMIP's and QUINOPT's internals
clear;
yalmip clear;
quinopt clear;

% Then,we set some problem parameters: the Reynolds number, the period Lambda_y
% in the y direction, the degree of the linear term PHI in the storage
% functional, and the maximum horizontal wavenumber to test:
Re       = 200;
Lambda_y = 4*pi;
PHI_DEG  = 15;
k_max    = 5;

% Finally, we define the integration variables, the flow variables, and the
% boundary conditions
z = indvar(0,1);
[u,w] = depvar(z);
BC = [u(0); u(1); w(0); w(1); w(0,1); w(1,1)];

%% Setting up the optimization variables
% The optimization variables are:
% 1) a, a scalar known in the literature as "balance parameter"
% 2) U, the upper bound on the time-average bulk dissipation
% 3) PHI, a polynomial of degree PHI_DEG satisfying the conditions PHI(0)=0 and
%    PHI(1)=0. We set up PHI as a polynomial in the Legendre basis, which is how
%    QUINOPT represents variable internally.
parameters a U
PHI = legpoly(z,PHI_DEG);

% The boundary conditions on PHI are specified in a vector, using the command
% legpoly() to compute the boundary values of PHI.
PHI_BC = [legpolyval(PHI,0)==0; legpolyval(PHI,1)==0];

% Finally, we will need the first two derivatives of PHI. These can be computed
% using the function jacobian().
D1PHI = jacobian(PHI,z);
D2PHI = jacobian(D1PHI,z);

%% Setting up the inequality constraints
% To set up the integral inequality constraints, we construct a vector EXPR 
% containing the integrand of each inequality.
n = 0;
k = 0;
EXPR(1) = (a-1)*u(z,1)^2 + D2PHI*u(z) + U-1;
while k<k_max
    n = n+1;
    k = 2*pi*n/Lambda_y;
    EXPR(end+1) = (a-1)*( u(z,1)^2 + k^2*u(z)^2 ) ...
                 +(a-1)*( w(z,2)^2/k^2 + 2*w(z,1)^2 + k^2*w(z)^2 ) ...
                 + Re*( a+D1PHI )*u(z)*w(z);
end

%% Solve the problem
% The aim of the optimization is to minimize the upper bound U on the energy
% dissipation. To specify the extra constraints on the optimization variables,
% we use the command quinopt() with five inputs: EXPR and BC to specify the
% integral inequalities, U as the objective, an options structure (empty) and
% the additional constraints PHI_BC.
quinopt(EXPR,BC,U,[],PHI_BC);
U_OPT = value(U);

% Print some solution info
fprintf('\nExample 9:\n')
fprintf('==============================\n')
fprintf('Reynolds number: %g\n',Re)
fprintf('Period Lambda_y: %g*pi\n',Lambda_y/pi)
fprintf('Optimal bound U: %g\n',U_OPT)
fprintf('==============================\n')

%% Inspecting the solution
% As a final step, we can inspect the solution to see what the optimal choice of
% the polynomial PHI is. To do so, we use the command plot(), which is
% overloaded on Legendre polynomials.
subplot(2,1,1)
plot(0:0.01:1,PHI)
xlabel('z'); ylabel('\phi(z)')

% In fact, to compare the optimization results obtained with QUINOPT to those in
% Plasting & Kerswell, J. Fluid Mech. 477, 363–379 (2003) it is convenient to
% use PHI and the balance parameter a to plot the usual "background field",
% given by a^{-1}*PHI(z) + z   
subplot(2,1,2)
plot(0:0.01:1,PHI/value(a)+z)
xlabel('z'); ylabel('a^{-1}\phi(z)+z')
drawnow;

%% END CODE