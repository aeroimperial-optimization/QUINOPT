%% example10.m
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

% Then we set some YALMIP options (cache the list of solvers, and run silently)
opts.YALMIP = sdpsettings('cachesolvers',1,'verbose',0); 

%% Set up the variables
% Then we set up the variables: the independent variable z in [0,1], the
% dependent variable T(z), and the Marangoni number M:
z = indvar(0,1);            
T = depvar(z);             
parameters M

% We also set up the boundary conditions, T(0)=0 and T'(1)=0:
BC = [T(0); T(1,1)];

%% Loop over wavenumbers
% To compute the maximum Marangoni number for stability at different wavenumbers
% k, we use a for loop. Before entering it, we define some useful variables:
% 1) k_val : the wavenumbers to consider in the loop
% 2) z_plot: coordinates for plotting
% 3) M_OPT : a vector to contain the optimal Marangoni numbers
% 4) Fk_deg: the polynomial degree to approximate the function Fk(z). Using
%            Fk_deg=10 is good enough accurate up to k=5.
k_val  = 1:0.25:5;                 
z_plot = 0:0.01:1;                 
M_OPT  = zeros(size(k_val));   
Fk_deg = 10; 

% Then we loop over each value of k, but first we display a header to print the
% results.
fprintf('\n=======================\n')
fprintf('|    k   |    M_OPT   |\n')
fprintf('=======================\n')
for i = 1:length(k_val)
    
    % First, extract the wavenumber for convenience:
    k = k_val(i);
    
    % The first step is to approximate the function Fk(z) with a polynomial
    % approximation, since QUINOPT only works on inequalities with polynomial
    % data. To do so, we define Fk(z) as a function handle, we compute its first
    % Fk_deg+1 Legendre expansion coefficients using the command "flt()" (Fast 
    % Legendre transform), and build a polynomial approximation FkPoly of degree
    %  Fk_deg using the command "legpoly()".
    Fk = @(z)k*sinh(k)/(sinh(2*k)-2*k).*(k*z.*cosh(k*z)-sinh(k*z)+(1-k*coth(k))*z.*sinh(k*z));
    leg_coef = flt(Fk,Fk_deg+1,[0,1]);            
    FkPoly = legpoly(z,Fk_deg,leg_coef);
    
    % We can easily plot the polynomial approximation vs the exact function
    % using the function "plot()", which is overloaded on polynomials built
    % using "legpoly()": 
    plot(z_plot,Fk(z_plot),'Linewidth',1.5); hold on
    plot(z_plot,FkPoly,'.','MarkerSize',12);
    xlabel('$z$','interpreter','latex','fontsize',12); 
    ylabel('$F_k(z)$','interpreter','latex','fontsize',12); 
    title(['F_k for k=',num2str(k)]);
    legend('Exact f_k(x)','Legendre series approximation','Location','southwest'); 
    axis([0 1 -0.2 0]); 
    drawnow; 
    hold off;
    
    % Having checked that the approximation is satisfactory, we can define the
    % integrand of the integral inequality and maximize M. 
    EXPR = T(z,1)^2 + k^2*T(z)^2 + M*FkPoly*T(z)*T(1);
    quinopt(EXPR,BC,-M,opts);
    M_OPT(i) = value(M);
    
    % Finally, we display the results
    fprintf('|  %4.2f  |  %8.4f  |\n',k,M_OPT(i))
    
end
fprintf('=======================\n\n')

%% Plot the results
clf; 
plot(k_val,M_OPT,'.-','linewidth',1,'markersize',12); 
xlabel('$k$','interpreter','latex','fontsize',12); 
ylabel('$M$','interpreter','latex','fontsize',12);

%% END CODE