clear; clc; close all
%% Set parameters 
L = 8;          % Lipschitz constant
M = 4;          % Smoothness constant
T = 300;        % Maximum iteration

mu = 0.01;      % Gradient estimation deviation upperbound
h = 0.1;        % Safety threshold
epsl = 1e-10;   % Convergence condition
rho = 0.9;      % Update rate of step length selection
c = 10^-4;      % Small constant in step length selection


options = struct;
options.L = L;
options.M = M;
options.mu = mu;
options.h = h;
options.epsl = epsl;
options.rho = rho;
options.c = c;

%% Define problem
x0 = [0,-4];                % Starting point
y_e = obj_fun(x0);          % Objective function
fi_e = fi_fun(x0);          % Constraint functions
x_hist = x0;                % Record x iteration
y_hist = y_e;               % Record y iteration
fi_hist = fi_e;             % Record fi iteration

gt_hist = [];               % Record gradient iteration

obj_hand = @obj_fun;        % Handle of objective function
fi_hand = @fi_fun;          % Handle of constraint function

%% Optimization loop
for iter = 1:T

    x_current = x_hist(end,:);

    [x_next,converged]=e_sls(x_current,obj_hand,fi_hand,options);

    y_next = obj_fun(x_next);
    fi_next = fi_fun(x_next);

    x_hist = [x_hist;x_next];
    y_hist = [y_hist;y_next];
    fi_hist = [fi_hist; fi_next];
  
    fprintf("Iter: %d | Obj: %4.2f | Max_Cons: %1.2f \n", iter,y_next,max(fi_next));
    
    if converged
        disp('converged');
        break;
    end

end

plot_figure(x_hist,y_hist,fi_hist)

%% Auxiliary functions

% Define objective function
function y = obj_fun(x)
    y = (x(1)-2.7)^2+0.5*(x(2)-0.5)^2-5;
end

% Define objective function
function fi = fi_fun(x)
    f1 = x(1)-2.7;
    f2 = -5-x(2);
    fi = [f1,f2];
end

% Plot figures
function []=plot_figure(x_hist,y_hist,fi_hist,O_fi,R_fi)
    n_itr = size(y_hist,1);
    
    figure(1)
    plot(linspace(0,n_itr,n_itr),y_hist(:,1),'b-',LineWidth=1); hold on
    line([0,size(y_hist,1)],[-5,-5],'Color','r','LineStyle','--','Linewidth',1);
    title('Objective');
    ylim([-6,20]);
    xlabel('$k$','Interpreter','latex')
    legend('objective')
    

    figure(2)
    plot(linspace(0,n_itr,n_itr),fi_hist(:,1),'b-'); hold on
    line([0,size(fi_hist,1)],[0,0],'Color','r','LineStyle','--','Linewidth',1);
    title('Constraints')
    legend('constraint')
    xlabel('$k$','Interpreter','latex')
     

    figure(3)
    plot(x_hist(:,1),x_hist(:,2),'b-',LineWidth=1); hold on    
    
    line([2.7,2.7],[-10,10],'Color','r','LineStyle','--','Linewidth',1);
    line([-10,10],[-5,-5],'Color','r','LineStyle','--','Linewidth',1);
    
    scatter(x_hist(1,1),x_hist(1,2),100,'ro','filled','MarkerFaceColor','r')
    scatter(2.7,0.5,100,'gp','filled','MarkerFaceColor','g')
    
    title('Optimization trajectory')
    xlabel('$x_1$','Interpreter','latex')
    ylabel('$x_2$','Interpreter','latex')
    legend('e-SLS','Constraint 1','Constraint 2','Start','Optimum','Location','northwest')
    xlim([-10,6])
    ylim([-8,6])
end