%% FeasibleSetApprox.m
% Example in Section 7.3
%
% Example to plot feasible region

%% Setup
clear; close all;
opts.YALMIP = sdpsettings('cachesolvers',1,'verbose',0);

dom = [-1 1];
bnds = [-5,5];

for N = 2
    
    fprintf('N = %i\n',N)
    figure();
    %hold on
    
    clearModel;
    yalmip clear
    
    x = indvar(dom(1),dom(2));
    [u,v] = depvar(x);
    parameters a b
    
    expr = u(x,1)^2 + v(x,1)^2 + a*(x^2)*u(x,1)*v(x,1) + 2*b*u(x)*v(x);
    BC = [u(dom(1)); u(dom(2)); v(dom(1)); v(dom(2))];
    
    
    %% Outer approximation
    % Use YALMIP's plot function, since easy problem
    fprintf('Finding outer approximation...')
    opts.rigorous = 0;
    cnstr = setQuadIntIneq(expr,BC,N,opts);
    plot([cnstr; bnds(1)<=a<=bnds(2); bnds(1)<=b<=bnds(2)],[a,b],'r')
    fprintf('done\n')
    
    %% Inner approximation
    % manual setup: lmi/plot screws up SOS constraints I think
    fprintf('Finding inner approximation...')
    hold on
    opts.rigorous = 1;
    opts.sosdeg = 16;
    [cnstr,data] = setQuadIntIneq(expr,BC,N,opts);
    
    % Compile SOS problem if needed (adapted from lmi/plot)
    if any(is(cnstr,'sos'))
        opts.YALMIP.sos.model = 2;
        cnstr = compilesos(cnstr,[],opts.YALMIP,data.SumOfSquaresParameters);
    end
    
    % Compile YALMIP problem (copied from yalmip/extras/export)
    [model,recoverdata,diagnostic,internalmodel] = export(cnstr,[],opts.YALMIP,[],[],0);
    if isempty(internalmodel) || (~isempty(diagnostic) && diagnostic.problem)
        error('Could not create model. Can you actually solve problems with this model?')
    end
    internalmodel.options.saveduals = 0;
    internalmodel.getsolvertime = 0;
    internalmodel.options.dimacs = 0;
    
    % Find local indices to set up cost (copied from YALMIP)
    pp = [a,b];
    localindex = [];
    for i = 1:2
        localindex = [localindex, find(ismember(recoverdata.used_variables,getvariables(pp(i))))];
    end
    
    x_opt = zeros(2,300);
    i = 1;
    for t = linspace(0,2*pi,300)
        c = [cos(t); sin(t)];
        internalmodel.c = 0*internalmodel.c;
        internalmodel.c(localindex) = c;
        sol  = feval(internalmodel.solver.call,internalmodel);
        xout = sol.Primal;
        x_opt(:,i) = xout(localindex(:));
        i = i+1;
    end
    fprintf('done\n\n')
    
    patch(x_opt(1,:),x_opt(2,:),'w')
    axis square; box on; axis([bnds(1),bnds(2),bnds(1),bnds(2)])
    set(gca,'FontSize',8);
    title(sprintf('$N = %i$',N),'Interpreter','latex','fontsize',8);
    xlabel('$\gamma_1$','Interpreter','latex','fontsize',10);
    ylabel('$\gamma_2$','Interpreter','latex','fontsize',10);
    
    set(gcf,'PaperUnits','centimeters')
    set(gcf,'PaperPosition',[0 0 3.75 3.75])
    figname = ['InnerOuterApprox_N',num2str(N),'.eps'];
    print(gcf,figname,'-depsc2','-r300')
    
end