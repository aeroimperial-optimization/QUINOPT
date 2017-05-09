%% example09.m
%
% In this example we compute the energy stability boundary of perturbation of 
% wavenumber k to the conduction state of Benard-Marangoni convection problem.
% This amounts to solving the optimization problem
%
%   maximize M
%               /1
%               |
%   subject to  | [ T'(z)^2 + k^2*T(z)^2 + M*Fk(z)*T(z)*T(1) ] dz >=0, T(0)=0=T'(1)
%               |
%               /0
%
% where
%
%         k sinh(k) (sinh(k z) + z sinh(k z) (k coth(k) - 1) - k z cosh(k z))
% Fk(z) = -------------------------------------------------------------------
%                           2 k - 2 cosh(k) sinh(k)
%
% We consider k = \pi for definiteness.
%
% This example demonstrates how to solve problems in which the unknown boundary 
% value of the dependent variables appear explicitly in the integral constraint.
% Moreover, it shows how QUINOPT can be used to approximate problems with data 
% that is non-polynomials using a truncated Legendre transform.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/04/2016
% Last Modified:    09/05/2017
% ----------------------------------------------------------------------- %

%% Initialization
% First, we clean up the workspace
clear;             
yalmip clear;      
quinopt clear;  


%% Set up the variables
% Then we set up the variables: the independent variable z in [0,1], the
% dependent variable T(z), and the Marangoni number M (which is the optimization
% variable).
z = indvar(0,1);            
T = depvar(z);             
parameters M

%% Construct a polynomial approximation to Fk(z)
% QUINOPT can only solve integral inequalities with polynomial data, but the 
% function :math:`F_k(z)` is clearly not a polynomial. Yet, QUINOPT can be used
% if we replace :math:`F_k(z)` with a polynomial approximation of degree d. To 
% do so, we first construct the exact function F_k(z) as a function handle, and 
% we compute its first d+1 Legendre expansion coefficients using the 
% fast Legendre transform command "flt()". Finally, we use these coefficients to
% build a polynomial approximation of degree d using the command "legpoly()". 
% Below, we set k=\pi and d=10 (this gives an accurate approximation at least up
% to k=5).
k = pi;
d = 10;
Fk = @(z)k*sinh(k)/(sinh(2*k)-2*k).*(k*z.*cosh(k*z)-sinh(k*z)+(1-k*coth(k))*z.*sinh(k*z));
leg_coef = flt(Fk,d+1,[0,1]);          % Compute the Legendre expansion coefficients
FkPoly = legpoly(z,d,leg_coef);        % Build a degree-d legendre polynomial with coefficients specified by leg_coef

% We can easily plot both Fk(z) and its polynomial approximation using the 
% function "plot()", which is overloaded on polynomials built using "legpoly()":
clf;                                                     % clear current figure
plot(0:0.01:1,Fk(0:0.01:1),'Linewidth',1.5); hold on;    % plot Fk(z)
plot(0:0.01:1,FkPoly,'.','MarkerSize',12); hold off;     % plot the polynomial approximation
xlabel('$z$','interpreter','latex','fontsize',12);
legend('F_k(z)','Polynomial approximation','Location','southwest');
axis([0 1 -0.2 0]);

%% Maximize M
% Once a polynomial approximation of Fk(z) has been constructed, the maximum 
% Marangoni number M satisfying the integral inequality at the top of the page 
% can be computed with QUINOPT. First, we define the integrand of the integral 
% inequality, and the boundary conditions on the dependent variable:
EXPR = T(z,1)^2 + k^2*T(z)^2 + M*FkPoly*T(z)*T(1);
BC = [T(0); T(1,1)];                                % The boundary conditions T(0)=0, T'(1)=0

%Then we maximize M by calling
quinopt(EXPR,BC,-M)
value(M)

% Note the negative sign in the objective function, which is needed because 
% QUINOPT minimizes the specified objective.

%% END CODE