%% example06.m
%
% Compute the maximum Grashoff number G such that a two-dimensional layer of
% fluid, driven by a surface shear stress, is stable using energy as the
% candidate Lyapunov function.
%
% This example introduces the use of multiple dependent variables.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/04/2016
% Last Modified:    28/04/2017
% ----------------------------------------------------------------------- %

%% Initialization
% First, we clean up
clear;          
yalmip clear;   
quinopt clear;

% Then, we set some problem variables: the horizontal period Lambda, the largest
% wavenumber to test.
lambda = 6*pi;
k_max  = 5; 

% To run in silent mode, we set YALMIP's option 'verbose' to 0. Also, we
% speed up the iterations by settin YALMIP's 'cachesolvers' option to 1.
% Finally, we set the degree of the Legendre expansion used internally by
% QUINOPT to N=10;
opts.YALMIP = sdpsettings('verbose',0,'cachesolvers',1);
opts.N = 15;


%% Loop over wavenumbers

% Display a header to print results
fprintf('\n================================\n')
fprintf('|   k    |    LB    |    UB    |\n')
fprintf('================================\n')

% Loop
k = 0;
n = 1;
while k <= k_max
    
        % Setup the variables
        y = indvar(0,1);   % define variable of integration over [0,1]
        [u,v] = depvar(y);  % define two dependent variables u and v
        parameters G        % define the optimization variable
        
        % Setup the inequality & the boundary conditions
        k(n) = 2*pi*n/lambda;
        expr = ( u(y,2)^2+v(y,2)^2 )/k(n)^2 ...
              +2*( u(y,1)^2+v(y,1)^2 ) ...
              +k(n)^2*( u(y)^2+v(y)^2 )  ...
              - G/k(n)*( u(y)*v(y,1) - u(y,1)*v(y) );
        bc = [u(0); u(1); u(0,1); u(1,2); ...        % BC on u
              v(0); v(1); v(0,1); v(1,2)];           % BC on v
        
        % Maximize G using an inner approximation of the integral inequality 
        % (we minimize -G, so we obtain a lower bound on the "true" optimal G)
        opts.method = 'inner';
        quinopt(expr,bc,-G,opts);
        LB(n) = value(G);
        
        % Maximize G using an outer approximation of the integral inequality 
        % (we minimize -G, so we obtain an upper bound on the "true" optimal G)
        opts.method = 'outer';
        quinopt(expr,bc,-G,opts);
        UB(n) = value(G);
        
        % Print progress
        fprintf('|  %4.2f  |  %6.2f  |  %6.2f  |\n',k(n),LB(n),UB(n))
        
        % Update variables for the next iteration. It is convenient to clear the
        % internal variables of QUINOPT and YALMIP to avoid accumulating unused 
        % internal variables that slows down the execution of the next iteration.
        quinopt clear;         % clear QUINOPT's internal variables
        yalmip clear;          % clear YALMIP's internal variables
        n = n+1;
        
end
fprintf('================================\n\n')

%% Plot the results
plot(k*lambda/2/pi,LB,'.-','displayname','lower bound on critical G'); hold on;
plot(k*lambda/2/pi,UB,'x-','displayname','upper bound on critical G'); hold off;
xlim([0, k_max*lambda/2/pi])
legend toggle
xlabel('k'); ylabel('UB and LB on optimal G');
%% END CODE