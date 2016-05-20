%% Example2.m
%
% Solve the Example in Section 8.2 of the reference paper.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/04/2016
% Last Modified:    12/05/2016
% ----------------------------------------------------------------------- %

%% PART 1: plot Figure

% Initialization
clear;          % clear workspace
xi = 1:0.1:10;  % values of xi to solve for
N = 10;         % Legendre truncation parameter

%To run in silent mode, we set YALMIP's option 'verbose' to 0. Also, we
%speed up the iteration by settin YALMIP's 'cachesolvers' option to 1.
opts.YALMIP = sdpsettings('verbose',0,'cachesolvers',1);

% Loop over values of xi to plot figure
gamma_cr1 = zeros(length(xi),1);
fprintf('PART 1: Looping over %i values of xi...',length(xi))
chrono = tic;
for i = 1:length(xi)

    % Setup variables
    x = indvar(-1,1);   % define variable of integration over [-1,1]
    [u,v] = depvar(x);  % define two dependent variables u and v
    parameters gamma    % define the optimization variable

    % Setup inequality
    expr = 16/xi(i)^2*( u(x,2)^2+v(x,2)^2 ) + 8*( u(x,1)^2+v(x,1)^2 ) + ...
        xi(i)^2*( u(x)^2+v(x)^2 )  - 2*gamma/xi(i)*( u(x)*v(x,1) - u(x,1)*v(x) );
    BC = [u(-1); u(1); u(-1,1); u(1,2)];        % BC on u
    BC = [BC; v(-1); v(1); v(-1,1); v(1,2)];    % BC on v

    % Solve optimization. To maximize gamma, we minimize -gamma (since
    % QUINOPT minimizes the objective function by default).
    % The Legendre series truncation parameter is fixed to N=10.
    quinopt(expr,BC,-gamma,[],[],10,opts);
    gamma_cr1(i) = value(gamma);

    % Clear internal variables of QUINOPT and YALMIP to avoid buildup of
    % unused variables that slows down the execution of the next iteration.
    clearModel;     % clear QUINOPT's internal variables
    yalmip clear;    % clear YALMIP's internal variables

end
fprintf('done in %.4f seconds\n',toc(chrono))

% Plot
plot(xi,gamma_cr1,'linewidth',1.5);
xlabel('$\xi$','interpreter','latex','fontsize',14);
% ylabel('$\gamma_\mathrm{cr}$','interpreter','latex','fontsize',14);

%% PART 2: check convergence with N for xi=pi
N = 4:4:20;
xi = [2:2:10,3.1469];

% Display a nice header
fprintf('PART 2: Increasing N for xi=pi\n')
fprintf('|============================================================================================================================|\n')
fprintf('|    |    xi = %5.2f     |    xi = %5.2f     |    xi = %5.2f     |    xi = %5.2f     |    xi = %5.2f     |    xi = %5.2f     |\n',xi)
fprintf('|============================================================================================================================|\n')
fprintf('|  N | gamma_cr   Time   | gamma_cr   Time   | gamma_cr   Time   | gamma_cr   Time   | gamma_cr   Time   | gamma_cr   Time   |\n')
fprintf('|============================================================================================================================|\n')
results = zeros(5,12);
for i=1:5
    for j=1:6
        
        % Setup variables
        x = indvar(-1,1);   % define variable of integration over [-1,1]
        [u,v] = depvar(x);  % define two dependent variables u and v
        parameters gamma    % define the optimization variable
        
        % Setup inequality
        expr = 16/xi(j)^2*( u(x,2)^2+v(x,2)^2 ) + 8*( u(x,1)^2+v(x,1)^2 ) + ...
            xi(j)^2*( u(x)^2+v(x)^2 )  - 2*gamma/xi(j)*( u(x)*v(x,1) - u(x,1)*v(x) );
        BC = [u(-1); u(1); u(-1,1); u(1,2)];        % BC on u
        BC = [BC; v(-1); v(1); v(-1,1); v(1,2)];    % BC on v
        
        % Solve optimization. To maximize gamma, we minimize -gamma (since
        % QUINOPT minimizes the objective function by default).
        chrono = tic;
        sol = quinopt(expr,BC,-gamma,[],[],N(i),opts);
        results(i,2*j) = toc(chrono);
        results(i,2*j-1) = value(gamma);
        
        % Clear internal variables of QUINOPT and YALMIP to avoid buildup of
        % unused variables that slows down the execution of the next iteration.
        clearModel;     % clear QUINOPT's internal variables
        yalmip clear;    % clear YALMIP's internal variables
       
   end 
end

fprintf('| %2.i | %8.4f %8.4f | %8.4f %8.4f | %8.4f %8.4f | %8.4f %8.4f | %8.4f %8.4f | %8.4f %8.4f |\n',[N',results]')
fprintf('|============================================================================================================================|\n')
