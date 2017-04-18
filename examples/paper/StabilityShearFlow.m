%% StabilityShearFlow.m
% Example in Section 7.1
%
% Motivating example: stability of a flow driven by a shear stress.
% See:  Motivating example, Section 3 of the paper.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/04/2016
% Last Modified:    17/04/2017
% ----------------------------------------------------------------------- %

%% PART 1: plot Figure

% Initialization
clear;          % clear workspace
quinopt clear;  % clear QUINOPT's internals
xi = 1:0.1:10;  % values of xi to solve for
savefigs = 0;   % 1 to save plots, 0 otherwise

%To run in silent mode, we set YALMIP's option 'verbose' to 0. Also, we
%speed up the iteration by settin YALMIP's 'cachesolvers' option to 1.
opts.YALMIP = sdpsettings('verbose',0,'cachesolvers',1);

for N = [3 6 9];  % Legendre truncation parameter
    
    % Loop over values of xi to plot figure
    gamma_cr_inn = zeros(length(xi),1);
    gamma_cr_out = zeros(length(xi),1);
    fprintf('Part 1, N = %i: Looping over %i values of xi...',N,length(xi))
    chrono = tic;
    
    for i = 1:length(xi)
        
        % Setup variables
        x = indvar(-1,1);   % define variable of integration over [-1,1]
        [u,v] = depvar(x);  % define two dependent variables u and v
        parameters gamma    % define the optimization variable
        
        % Setup inequality
        expr = 16/(xi(i)^2)*( u(x,2)^2+v(x,2)^2 ) + 8*( u(x,1)^2+v(x,1)^2 ) + ...
            xi(i)^2*( u(x)^2+v(x)^2 )  - 2*gamma/xi(i)*( u(x)*v(x,1) - u(x,1)*v(x) );
        BC = [u(-1); u(1); u(-1,1); u(1,2)];        % BC on u
        BC = [BC; v(-1); v(1); v(-1,1); v(1,2)];    % BC on v
        
        % Solve optimization. To maximize gamma, we minimize -gamma (since
        % QUINOPT minimizes the objective function by default).
        % The Legendre series truncation parameter is fixed to N=12.
        
        % Inner approximation (maximum lower bound on optimal gamma)
        opts.method = 'inner';
        opts.N = N;
        sol = quinopt(expr,BC,-gamma,opts);
        if sol.FeasCode==1
            % Infeasible for gamma~=0
            gamma_cr_inn(i) = 0;
        elseif sol.FeasCode==2
            % Unbounded
            gamma_cr_inn(i) = +Inf;
        else
            gamma_cr_inn(i) = value(gamma);
        end
        
        % Outer approximation (minimum upper bound on optimal gamma)
        opts.method = 'outer';
        sol = quinopt(expr,BC,-gamma,opts);
        if sol.FeasCode==1
            % Infeasible for gamma~=0
            gamma_cr_out(i) = 0;
        elseif sol.FeasCode==2
            % Unbounded
            gamma_cr_out(i) = +Inf;
        else
            gamma_cr_out(i) = value(gamma);
        end
        
        % Clear internal variables of QUINOPT and YALMIP to avoid buildup of
        % unused variables that slows down the execution of the next iteration.
        quinopt clear;         % clear QUINOPT's internal variables
        yalmip clear;          % clear YALMIP's internal variables
          
    end
    fprintf('done in %.4f seconds\n',toc(chrono))
    
    % Plot
    figure();
    plot(xi,gamma_cr_inn,'-','linewidth',1.5,'displayname','Inner Approximation'); hold on;
    plot(xi,gamma_cr_out,'-.','linewidth',1.5,'displayname','Outer Approximation'); box on;
    xlabel('$\xi$','interpreter','latex','fontsize',10);
    title(sprintf('$N = %i$',N),'interpreter','latex')
    axis([1,10,0,400])
    set(gca,'fontsize',8); ax = gca; ax.YAxis.Exponent = 2;
    
    
    % Save figures
    if savefigs
        set(gcf,'PaperUnits','centimeters')
        set(gcf,'PaperPosition',[0 0 5.0 4.5])
        figname = ['ShearFlowExample_N',num2str(N),'.eps'];
        print(gcf,figname,'-depsc2','-r300')
    end
    
end

%% PART 2: check convergence with N for xi=pi
N = 3:3:12;
xi = [3 6 9];

% Display a nice header
fprintf('PART 2: Increasing N for xi=pi\n')
fprintf('|================================================================================================================|\n')
fprintf('|    |            xi = %5.2f             |             xi = %5.2f            |             xi = %5.2f            |\n',xi)
fprintf('|================================================================================================================|\n')
fprintf('|  N |   Inner    Time    Outer    Time  |   Inner    Time    Outer    Time  |   Inner    Time    Outer    Time  |\n')
fprintf('|================================================================================================================|\n')
results = zeros(length(N),12);
for i=1:length(N)
    for j=1:length(xi)
        
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
        
        % Inner approximation
        opts.method = 'inner';
        chrono = tic;
        opts.N = N(i);
        sol = quinopt(expr,BC,-gamma,opts);
        results(i,4*j-2) = toc(chrono);
        if sol.FeasCode==1
            % Infeasible for gamma~=0
            results(i,4*j-3) = 0;
        elseif sol.FeasCode==2
            % Unbounded
            results(i,4*j-3) = Inf;
        else
            results(i,4*j-3) = value(gamma);
        end
        
        % Outer approximation
        opts.method = 'outer';
        chrono = tic;
        sol = quinopt(expr,BC,-gamma,opts);
        results(i,4*j) = toc(chrono);
        if sol.FeasCode==1
            % Infeasible for gamma~=0
            results(i,4*j-1) = 0;
        elseif sol.FeasCode==2
            % Unbounded
            results(i,4*j-1) = Inf;
        else
            results(i,4*j-1) = value(gamma);
        end
        
        
        % Clear internal variables of QUINOPT and YALMIP to avoid buildup of
        % unused variables that slows down the execution of the next iteration.
        quinopt clear;   % clear QUINOPT's internal variables
        yalmip clear;    % clear YALMIP's internal variables
        
    end
end

fprintf('| %2.i | %8.4f %6.2f %10.4f %6.2f | %8.4f %6.2f %10.4f %6.2f | %8.4f %6.2f %10.4f %6.2f |\n',[N',results]')
fprintf('|================================================================================================================|\n')
