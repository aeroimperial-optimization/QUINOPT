%% AdditionalEx5.m
%
% Compute the energy stability boundary of perturbation of wavenumber k to
% conduction state of Benard-Marangoni convection problem. This example
% demonstrates how to solve problems in which the unknown boundary value of
% the dependent variables appear explicitly in the integral constraint.
%
% NOTE: Strictly speaking, this problem involves non-polynomial coefficients.
% We approximate these with a truncated Legendre series (we checked that 
% enough terms are used).

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    08/04/2016
% Last Modified:    12/05/2016
% ----------------------------------------------------------------------- %

%% CODE

clear; clearModel;         
figure(gcf); clf;
time = tic;

%% Variables
x = indvar(0,1);            % independent variable in [0,1]
T = depvar(x);              % temperature fluctation
parameters Ma               % Variable Marangoni number
BC = [T(0); T(1,1)];        % BC

%% Loop over wavenumber
approxDeg = 20;                 % 20 should be accurate up to k=20
kval = 1:0.1:5;                 % values of k to test
xval = 0:0.01:1;                % x coordinate for plotting
Ma_cr = zeros(length(kval),1);  % vector of critical Ma

% Set options to make YALMIP faster - see YALMIP's documentation!
opts.YALMIP = sdpsettings('cachesolvers',1,'verbose',0);   

for i = 1:length(kval)
    
    % Define function fk, find its Legendre series approximation and plot
    k = kval(i);
    fk = @(z)k*sinh(k)/(2*(sinh(k)*cosh(k)-k)).*(k*z.*cosh(k*z)-sinh(k*z)+(1-k*coth(k))*z.*sinh(k*z));
    fcoef = flt(fk,approxDeg+1,[0,1]);            
    F = legpoly(x,approxDeg,fcoef);               % define function F as a legpoly
    figure(gcf); plot(xval,fk(xval),'Linewidth',1.5); hold on
    plot(xval,legpolyval(F,xval),'.','MarkerSize',12)
    xlabel('x'); ylabel('f_k(x)'); title(['f_k for k=',num2str(k)]);
    legend('Exact f_k(x)','Legendre series approximation','Location','southwest'); 
    axis([0 1 -0.2 0]); drawnow; hold off;
    
    % Set quadratic form
    Q = T(x,1)^2 + k^2*T(x)^2 + Ma*F*T(x)*T(1);
    
    % Solve
    quinopt(Q,BC,-Ma,[],[],[],opts);
    Ma_cr(i) = value(Ma);
    
end
time = toc(time);

% display & plot
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++')
fprintf('QUINOPT example 5: Solved %i wave numbers in %4.2f seconds.\n',length(kval),time);
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++'); disp(' ')
clf; plot(kval,Ma_cr,'.-','linewidth',1,'markersize',12); 
xlabel('k','fontsize',12); ylabel('Ma','fontsize',12)