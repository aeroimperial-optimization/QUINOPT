%% Example7.m
%
% Find values a, b, c, d to minimize a+b+c+d such that
%
% /1
% |
% |  ( u_x^2 + a*u*u_x + u^2 + b*u + c*u_x + d )dx  >= 0
% |
% /0
%
% holds for all functions u(x) subject to the Dirichlet boundary conditions 
% u(0)=0.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    17/04/2017
% Last Modified:    17/04/2017
% ----------------------------------------------------------------------- %


%% Initialization
%
% First, we remove any existing variables to avoid unexpected dependencies
clear;          % clear workspace
yalmip clear;   % clear YALMIP's internal variables
quinopt clear;  % clear QUINOPT's internal variables

% Then we initialize the independent variable of integration, the dependent
% variable u and the optimization variable gamma with the commands
x = indvar(0,1);        % define independent variable with domain [-1,1]
u = depvar(x);          % define u, which depends on x
parameters a b c d e    % define the optimsation variables

%% Construct the problem
% The integrand of the integral inequality is constructed with the commands
expr = u(x,1)^2 + a*u(x)*u(x,1) + u(x)^2 + b*u(x) + c*x^2*u(x,1) + d*u(0) + e + d ;
%
% where u(x,1) specifies the first derivative of u evaluated at the
% variable point x. 

% The boundary conditions u(-1)=0 and u(1)=0 are specified by constructing 
% the vector
BC = [u(0)];
%
% where the syntax u(0) specifies the value of $u$ at the boundary point 
% x = -1.

% Finally we specify the objective to minimize as
obj = a + b + c + d + e;

%% Solve the problem 
% We solve the problem for a number of values of theLegendre series 
% truncation parameter N. 

% To run in silent mode, we set YALMIP's option 'verbose' to 0. Also, we
% speed up the iteration by settin YALMIP's 'cachesolvers' option to 1.
opts.YALMIP = sdpsettings('verbose',0,'cachesolvers',1,'solver','sdpt3');
opts.rigorous = 0;

for N = 1:10
    
    % Call QUINOPT
    quinopt(expr,BC,obj,[],[],N,opts);
    fprintf('[a b c d e] = [%f %f %f %f %f]',value(a),value(b),value(c),value(d),value(e))
    fprintf('  obj = %f\n',value(obj))

end

%% END SCRIPT
