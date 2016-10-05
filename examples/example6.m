%% AdditionalEx6.m
% 
% We study the weighted L2 stability of a linear PDE
%
% u_t = u_xx+ (k - 24*x + 24*x^2)*u,  
%
% with boundary conditions
%
% u(-1)=u(1)=0
%
% The goal is to maximize k. The maximum k should be given by energy 
% stability (linear stability analysis gives the same eigenvalue problem 
% obtained by studying the energy stability using the calculus of variations).
% This example shows how QuadIntIneq automatically performs integration by
% parts when needed to obtain a sensible problem.
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
% Last Modified:    12/05/2016
% ----------------------------------------------------------------------- %

%% Clear previous model and variables
clearModel
clear

%% Problem variables
x = indvar(0,1);
u = depvar(x);
bc = [u(0); u(1)];

%% Classic energy stability: find the optimal k
% Check if 
%
%       /1
% V(u)= |  0.5*p(x)*u^2 dx 
%       /-1 
%
% is a Lyapunov function by testing if 
%
%       /1
% V_t(u)= |  p(x)*u*u_t dx  <= 0 for all u.
%       /-1 
%
% Use automatic integration by parts to formulate a well-posed inequality.
sdpvar k
u_t = u(x,2) + (k-24*x+24*x^2)*u(x);      % the right-hand side of the PDE
expr = u(x)*u_t;                          % the integrand

options.YALMIP = sdpsettings('verbose',0);
quinopt(-expr,bc,-k,[],[],5,options);   % require positivity of -V_t(u)!
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('QUINOPT example 6')
fprintf('Optimal k from energy stability: k = %4.4f\n',value(k))

%% Weighted energy stability
% Test if a weighted energy Lyapunov function achieves the same energy
% stability limit: this is not the case when the SOS method of [1] is used.
% Specifically, we use 
%
%       /1
% V(u)= |  0.5*p(x)*u^2 dx 
%       /-1
%
% as a Lyapunov function, where p(x) is a polynomial of degree d. 
% For a constant c>0 we therefore require:
%
% 1) V(u) >= c*||u||^2
% 2) V_t(u) = \int_{-1}^{1} p(x)*u*u_t dx <= 0.
%
% Since both expressions are homogeneous in p(x), take c=1 without loss of
% generality. Moreover, we require that the sum of the coefficients of the
% monomials x, x^2, ..., x^d in p(x) is larger than 1, so we don't get a 
% constant polynomial p(x) (this corresponds to the usual energy stability 
% choice).

% Define p using YALMIP's "polynomial" function
d = 10;
[p,c] = polynomial(x,d);      % the Lyapunov function is \int p(x)*u^2 dx
cnstr = sum(c(2:end))>=1;     % avoid solution p(x)=const.

% Setup feasibility problem for optimal value of k (actually, we try 0.9999
% of it to account for roundoff errors that might make the optimal k
% obtained from the energy stability slightly infeasible)
k = 0.9999*value(k);        
u_t = u(x,2) + (k-24*x+24*x^2)*u(x);
expr = [ p*u(x)^2; ...                    % positivity of V(u)
        -p*u(x)*u_t];                     % positivity of -V_t(u)

% Solve the feasibility problem. We need to specify that the coefficients c
% and d are parameters and not independent variables to obtain the correct
% answer.
sol = quinopt(expr,bc,[],cnstr,[],[],options);

% Plot p(x) if problem was successfully solved
if sol.FeasCode==0
    ttot = sol.solutionTime+sol.setupTime;
    fprintf('Feasible nonconstant p(x) could be found in %4.2f sec.\n',ttot);
    cval = clean(value(c),1e-8);        % YALMIP: remove small coefficients
    xp = 0:0.01:1;
    pval = polyval(flipud(cval),xp);    % polyval uses reverse order of coefficients
    figure(gcf); clf; plot(xp,pval,'Linewidth',1.5); xlabel('x'); ylabel('p(x)')
else
    fprintf('Feasible nonconstant p(x) could not be found (FeasCode = %2i).\n',sol.FeasCode)
end
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++'); disp(' ')

