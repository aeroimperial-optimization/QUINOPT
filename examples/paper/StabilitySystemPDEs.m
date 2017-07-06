%% StabilitySystemPDEs.m
% Example in Section 7.2
%
% Study linear, energy and weighted L2 stability of a linear PDE system
%
% u_t = gamma*u_xx+ u + 1.5*v,
% v_t = gamma*v_xx+ 5*u + 0.2*v,
%
% with Dirichlet boundary conditions
%
% u(0)=u(1)=v(0)=v(1)=0
%
% The goal is to minimize gamma (with a line search to avoid nonconvexity)
% such that the functional
%
%           /1
% V(u,v) =  |   [u v] [P11(x) P12(x)] [u]  dx
%           |         [P12(x) P22(x)] [v]
%           /0

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    11/04/2016
% Last Modified:    17/04/2017
% ----------------------------------------------------------------------- %

%% Initialization
clear; 
quinopt clear;
a=1; b=1.5; c=5; d=0.2;  % problem variables

%% Linear stability - find critical value with bisection method
% Set Chebyshev differentiation matrices
N = 64;
xcheb = cos(pi*(0:N)/N)';
r = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(xcheb,1,N+1);
dX = X-X';
D = (r*(1./r)')./(dX+eye(N+1));
D = D-diag(sum(D,2));
D2 = D^2;
D2 = D2(2:end-1,2:end-1);   % enforce BC
I = eye(N-1);

% Solve linear eigenvalue problem (use good starting guess from preliminary experiments)
gamma1 = 0.33;
gamma2 = 0.4;
gamma_lin = 0.35;
ecr = 1;
while abs(ecr)>1e-10
    A = [4*gamma_lin*D2+a*I, b*I; c*I, 4*gamma_lin*D2+d*I];
    e = sort(eig(-A));
    ecr = e(1);
    if ecr > 0
        gamma2 = gamma_lin;
    else
        gamma1 = gamma_lin;
    end
    gamma_lin = (gamma1+gamma2)/2;
end

% Display
fprintf('Linear stability  : gamma_lin = %8.6f\n',gamma_lin);

%% Energy stability
% This problem is convex so we can compute gamma_cr with optimization.

% Display nice header
fprintf('Lyapunov stability:\n')
fprintf('|===================================================|\n')
fprintf('| deg(P) |   gamma_cr   |  time (s)  | Av. it. time |\n')
fprintf('|===================================================|\n')

% Problem variables
x = indvar(0,1);
[u,v] = depvar(x);
parameters gamma
U = [u(x); v(x)];                         % the state vector
Ut = [gamma*u(x,2) + a*u(x) + b*v(x);...
    gamma*v(x,2) + c*u(x) + d*v(x)];    % the right-hand side of the PDE

% Set and solve inequality with N=10 (gives converged answer). To run in
% silent mode, we set YALMIP's option 'verbose' to 0. Also, we speed up
% YALMIP by setting the 'cachesolvers' option to 1.
P = eye(2);                         % P = identity matrix
expr = -U'*P*Ut;                    % time derivative of V is negative
bc = [u(0); u(1); v(0); v(1)];      % boundary conditions
opts.YALMIP = sdpsettings('verbose',0,'cachesolvers',1);

% Inner approximation
time = tic;
opts.N = 10;
quinopt(expr,bc,gamma,opts);       % require positivity of -V_t(u)!
time = toc(time);

% Display
gamma_en = value(gamma);
fprintf('| 0(P=I) |   %8.6f   |  %8.4f  |    ------    | (inner approx.)\n',gamma_en,time)

% Outer approximation
opts.method = 'outer';
time = tic;
quinopt(expr,bc,gamma,opts);       % require positivity of -V_t(u)!
time = toc(time);

% Display
gamma_en = value(gamma);
fprintf('| 0(P=I) |   %8.6f   |  %8.4f  |    ------    | (outer approx.)\n',gamma_en,time)

%% Weighted energy stability
% Do line search with bisection method to avoid nonconvexity. Try to locate
% the value of gamma_cr at which the problem becomes infeasible.
for degp = 0:2:6;
    
    % ------------------------------------------------------------------------ %
    for rig = [1,0]
        % 1 = inner approx.
        % 0 = outer approx.
        if rig
            opts.method = 'inner';
        else
            opts.method = 'outer';
        end
        
        gamma1 = 0.4;                   % surely feasible (larger than gamma_en)
        gamma2 = 0.34;                  % surely infeasible (smaller than gamma_lin)
        gamma_cr = gamma_lin;           % start from gamma_lin!
        
        time = tic;
        itAverageTime = 0;
        while gamma1-gamma2 > 1e-3
            % Clear internal variables of QUINOPT and YALMIP to avoid buildup of
            % unused variables that slows down the execution of the next iteration.
            yalmip clear;      % clear YALMIP's internal variables
            quinopt clear;     % clear QUINOPT's internal variables
            
            % Problem variables
            x = indvar(0,1);
            [u,v] = depvar(x);
            parameters gamma
            
            % Define P as symmetric polynomial matrix
            P = polyMat(x,degp,[2 2],'symm');
            
            % Setup feasibility problem
            U = [u(x);v(x)];                               % the state vector
            Ut = [gamma_cr*u(x,2) + a*u(x) + b*v(x);...
                gamma_cr*v(x,2) + c*u(x) + d*v(x)];         % the right-hand side of the PDE
            expr = [ U'*(P-eye(2))*U; ...                  % positivity of V(u)
                -U'*P*Ut];                             % positivity of -V_t(u)
            bc = [u(0); u(1); v(0); v(1)];                 % boundary conditions
            
            % Solve the feasibility problem with N=10
            opts.N = 10;
            chron = tic;
            sol = quinopt(expr,bc,[],opts);
            itAverageTime = (itAverageTime + toc(chron))/2;
            
            % Check if problem was successfully solved and update gamma_cr
            if sol.FeasCode==0 || sol.FeasCode==-2
                gamma1 = gamma_cr;
            else
                gamma2 = gamma_cr;
            end
            gamma_cr = 0.5*(gamma2+gamma1);
        end
        time = toc(time);
        
        % Display
        if rig == 1
            fprintf('|   %i    |   %8.6f   |  %8.4f  |    %6.4f    | (inner approx.)\n',degp,gamma1,time,itAverageTime)
        else
            fprintf('|   %i    |   %8.6f   |  %8.4f  |    %6.4f    | (outer approx.)\n',degp,gamma1,time,itAverageTime) 
        end
        
    end
end
fprintf('|===================================================|\n')
