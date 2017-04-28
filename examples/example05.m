%% example05.m
% 
% We study the nonlinear stability of a linear PDE
%
% u_t = u_xx+ (k - 24*x + 24*x^2)*u,  
%
% with boundary conditions
%
% u(-1)=u(1)=0
%
% by looking for a Lyapunov function of the form
%    
%        /1
%        |
% V(u) = | p(x)*u(x)^2 dx
%        |
%        /0
%
% The polynomial p(x) should be chosen such  that V is positive definite, and
% such that the time derivative of V is negative semidefinite.
%
% This example was taken from: 
% [1] Valmorbida, Ahmadi, Papachristodoulou, "Semi-definite programming and 
% functional inequalities for Distributed Parameter Systems", 53rd IEEE 
% Conference on Decision and Control, 2014.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/04/2016
% Last Modified:    28/04/2017
% ----------------------------------------------------------------------- %

%% Clear previous model and variables
clear;             % clean workspace 
yalmip clear;      % clear YALMIP's internal variables
quinopt clear;     % clear QUINOPT's internal variables

%% Define the problem variables
% We define the integration variable, the dependent variable, its boundary
% condition
x = indvar(0,1);
u = depvar(x);
bc = [u(0); u(1)];

% We then define the right-hand side of the PDE:
k = 15;
u_t = u(x,2) + (k-24*x+24*x^2)*u(x);

% Finally, we define the polynomial kernel of the Lyapunov function. We use a
% polynomial in the Legendre basis, which is how QUINOPT represents the
% variables internally.
p = legpoly(x,4);      

%% Define Lyapunov's inequalities
% We want to check if
%
%       /1
% V(u)= |  0.5*p(x)*u^2 dx 
%       /-1
%
% as a Lyapunov function. Lyapunov's theorem tells us that the PDE is stable
% around the zero zolution if there exists a constant c>0 such that
%
% 1) V(u) >= c*||u||^2
% 2) V_t(u) = \int_{-1}^{1} p(x)*u*u_t dx <= 0.
%
% Since both expressions are homogeneous in p(x), we may take c=1 without loss 
% of generality.
%
% To set up the inequalities in QUINOPT, we specify their integrands as the
% elements of a vector EXPR
EXPR = [ (p-1)*u(x)^2; ...                % V(u) - ||u||^2 >= 0
        -p*u(x)*u_t];                     % -V_t(u) >=0.
%
% Note the minus sign in the second integrand.
    
%% Solve for p(x)
% To solve for p(x), we use the command quinopt() with one output argument:
SOL = quinopt(EXPR,bc);
%
% The output SOL contains some information about the solution. In particular,
% the field SOL.FeasCode helps veryfying that the problem was successfully
% solved (in this case, SOL.FeasCode==0). Run 
%
%   >> help quinopt
%   >> quinoptFeasCode 
%
% for more information on QUINOPT's outputs, and the meaning of feasibility 
% codes (these are the same ones used by YALMIP).
%
% Try changing the value of k: for which value does a suitable p(x) stop
% existing?

%% Plot and display p(x) if successful
% Plot a polynomial defined using the command "legpoly()" is simple, because the
% usual plotting function "plot()" is overloaded. To display it, we use YALMIP's
% commands "value()" and "sdisplay()", which are overloaded on polynomials
% defined with "legpoly()".
fprintf('\n=================================================================\n')
fprintf('Solution report:\n')
fprintf('----------------\n')
if SOL.FeasCode==0
    ttot = SOL.solutionTime+SOL.setupTime;
    fprintf('A feasible p(x) was found in %4.2f seconds:\n',ttot);
    s = sdisplay(value(p));
    fprintf('p(x) = %s\n',s{1});
    clf
    plot(0:0.01:1,p)
    xlabel('x'); ylabel('p(x)');
else
    fprintf('A feasible p(x) could not be found (FeasCode = %2i).\n',SOL.FeasCode)
end
fprintf('=================================================================\n\n')

