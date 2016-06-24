%% StabilityShearFlow.m
% Example in Section 7.1
%
% Motivating example: stability of a flow driven by a shear stress.
% See:  Motivating example, Section 3 of the paper.
%       Example in Section 8.1 of the paper.
%
% NOTE: For N<4, the outer approximation displays errors - this is due to the
%       fact that for N<4 there are no polynomial satisfying the boundary
%       conditions except for the zero polynomial, and therefore the outer
%       SDP relaxation becomes an unconstrained minimization. The error messages
%       should be interpreted as an indication that the problem is unbounded
%       from below. The exact error message depends on the solver used. Examples
%       are:
%
%       MOSEK:
%         *** Error(1201): prob.blc has invalid dimension
%
%       SeDuMi:
%           SeDuMi had unexplained problems, maybe due to linear dependence?
%           YALMIP tweaks the problem (adds 1e6 magnitude bounds on all variables) and restarts...

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
savefigs = 0;   % 1 to save plots, 0 otherwise

%To run in silent mode, we set YALMIP's option 'verbose' to 0. Also, we
%speed up the iteration by settin YALMIP's 'cachesolvers' option to 1.
opts.YALMIP = sdpsettings('verbose',0,'cachesolvers',1);

for N = [4 8 12];  % Legendre truncation parameter
    
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
        expr = 16/xi(i)^2*( u(x,2)^2+v(x,2)^2 ) + 8*( u(x,1)^2+v(x,1)^2 ) + ...
            xi(i)^2*( u(x)^2+v(x)^2 )  - 2*gamma/xi(i)*( u(x)*v(x,1) - u(x,1)*v(x) );
        BC = [u(-1); u(1); u(-1,1); u(1,2)];        % BC on u
        BC = [BC; v(-1); v(1); v(-1,1); v(1,2)];    % BC on v
        
        % Solve optimization. To maximize gamma, we minimize -gamma (since
        % QUINOPT minimizes the objective function by default).
        % The Legendre series truncation parameter is fixed to N=12.
        
        % Inner approximation (upper bound)
        opts.rigorous = 1;
        quinopt(expr,BC,-gamma,[],[],N,opts);
        gamma_cr_inn(i) = value(gamma);
        
        % Inner approximation (lower bound)
        opts.rigorous = 0;
        quinopt(expr,BC,-gamma,[],[],N,opts);
        gamma_cr_out(i) = value(gamma);
        
        % Clear internal variables of QUINOPT and YALMIP to avoid buildup of
        % unused variables that slows down the execution of the next iteration.
        clearModel;         % clear QUINOPT's internal variables
        yalmip clear;       % clear YALMIP's internal variables
        
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
N = 4:4:40;
xi = [3 6 9];

% Display a nice header
fprintf('PART 2: Increasing N for xi=pi\n')
fprintf('|===================================================================================================================|\n')
fprintf('|    |            xi = %5.2f              |             xi = %5.2f             |             xi = %5.2f             |\n',xi)
fprintf('|===================================================================================================================|\n')
fprintf('|  N |   Inner    Time     Outer    Time  |   Inner    Time     Outer    Time  |   Inner    Time     Outer    Time  |\n')
fprintf('|===================================================================================================================|\n')
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
        opts.rigorous = 1;
        chrono = tic;
        sol = quinopt(expr,BC,-gamma,[],[],N(i),opts);
        results(i,4*j-2) = toc(chrono);
        results(i,4*j-3) = value(gamma);
        
        % Outer approximation
        opts.rigorous = 0;
        chrono = tic;
        sol = quinopt(expr,BC,-gamma,[],[],N(i),opts);
        results(i,4*j) = toc(chrono);
        results(i,4*j-1) = value(gamma);
        
        % Clear internal variables of QUINOPT and YALMIP to avoid buildup of
        % unused variables that slows down the execution of the next iteration.
        clearModel;     % clear QUINOPT's internal variables
        yalmip clear;    % clear YALMIP's internal variables
        
    end
end

% fprintf('| %2.i | %8.4f %6.2f %10.4f %6.2f | %8.4f %6.2f %10.4f %6.2f | %8.4f %6.2f %10.4f %6.2f  |\n',[N',results]')
fprintf(' %2.i & %8.4f & %6.2f & %10.4f & %6.2f & %8.4f & %6.2f & %10.4f & %6.2f & %8.4f & %6.2f & %10.4f & %6.2f  \\\\\n',[N',results]')
fprintf('|===================================================================================================================|\n')
