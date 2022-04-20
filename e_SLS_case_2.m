clear; clc; close all

%% Define e-SLS parameters

L = 8;                      % Lipschitz constant
M = 4;                      % Smoothness constant
T = 300;                    % Maximum iteration

mu = 0.01;                  % Gradient estimation deviation upperbound
h = 0.05;                   % Safety threshold
epsl = 1e-10;               % Convergence condition
rho = 0.9;                  % Update rate of step length selection
c = 10^-4;                  % Small constant in step length selection
use_newton_direction = 1;   % Use quasi-newton direction or steepest descent: 
                            % 1-quasi_newton, 0-steepest descent

%% Define problem
x0 = [-2,-1];               % Starting point
y_0 = obj_fun(x0);          % Objective function
fi_0 = fi_fun(x0);          % Constraint functions
x_hist = x0;                % Record x iteration
y_hist = y_0;               % Record y iteration
fi_hist = fi_0;             % Record fi iteration

gt_hist = [];               % Record gradient iteration


d = size(x0,2);             % Dimension of the problem
m = size(fi_0,2);           % Number of constraints

H = eye(d);

%% Optimization loop
for iter = 1:T

    x_current = x_hist(end,:);
    y_current = y_hist(end,:);
    fi_current = fi_hist(end,:);
    lambda_sol = [];

    % Compute step length for finite difference
    vk = min(2*mu/(sqrt(d)*M), min(-fi_current)/(2*L));   

    X_current = repmat(x_current,d,1);
    VT = eye(d)*vk;
    X_current = X_current+VT;
    Y_current = [];
    FI_current = [];
    for i = 1:d
        y_e = obj_fun(X_current(i,:));
        fi_e = fi_fun(X_current(i,:));
        Y_current = [Y_current; y_e];
        FI_current = [FI_current; fi_e];
    end
    
    G0 = (Y_current - y_current)/vk;        % Compute gradient estimator for objective
    GI = (FI_current - fi_current)/vk;      % Compute gradient estimator for constraints
    gt_hist = [gt_hist; G0'];  

    if use_newton_direction
        if iter>1
            x_last = x_hist(end-1,:);
            gt_last = gt_hist(end-1,:)';
            x_diff = x_current'-x_last';
            gt_diff = G0 - gt_last;
            
            psi = 1/(gt_diff'*x_diff);

            H = (eye(d) - psi*x_diff*gt_diff')*H*(eye(d)-psi*gt_diff*x_diff') + psi*(x_diff*x_diff');
        end

        p = -H*G0;
    else
        p = -G0;                                % Define the search direction
    end

    % Gradient estimation deviation
    Delta_ub = sqrt(d)*vk*M/2;
    
    if norm(p)~=0
        [~,I]=max(p~=0);                        % The first nonzero element
        e = eye(d);
        e = e(:,I);
        GI_mod = Delta_ub*norm(p)*e/(e'*p);     % Modification on GI
        GI_hat = GI + GI_mod;
    else
        GI_hat = GI;
    end

    % Search direction projection if min(-fi) <= safety threshold
    if min(-fi_current) <= h
        p_orig = p;

        % Find active constraint indices
        A = find(-fi_current<=h);
        % Solve non-negative least squares problem
        [lambda_sol,resnorm,resvec] = lsqnonneg(GI_hat(:,A),p_orig);    

        p_proj = -resvec;
        p = p_proj;

        % Recompute GI hat
        if norm(p)~=0
            [~,I]=max(p~=0);                        % The first nonzero element
            e = eye(d);
            e = e(:,I);
            GI_mod = Delta_ub*norm(p)*e/(e'*p);     % Modification on GI
            GI_hat = GI + GI_mod;
        else
            GI_hat = GI;
        end
    end

    % Centers and radii of each safe set of fi
    O_fi = x_current-GI_hat'./M;            
    R_fi = sqrt((vecnorm(GI_hat',2,2)/M).^2-2*fi_current'/M);

    % Compute the upper bound of the step length
    alpha_hat = min(R_fi);
    while min(R_fi - vecnorm((x_current + alpha_hat*p') - O_fi,2,2)) <= 0
        alpha_hat = alpha_hat*rho;
    end
    alpha = alpha_hat;

    y_k = y_current;
    fi_k = fi_current;

    y_k1 = obj_fun(x_current+alpha*p');
    fi_k1 = fi_fun(x_current+alpha*p');

    % Select safe step length
    while y_k1 > y_k + c*alpha*(G0'*p)
        alpha = alpha* rho;

        y_k1 = obj_fun(x_current+alpha*p');
        fi_k1 = fi_fun(x_current+alpha*p');
    end

    % Update iteration
    
    x_next = x_current + alpha*p';
    y_next = obj_fun(x_next);
    fi_next = fi_fun(x_next);

    x_hist = [x_hist;x_next];
    y_hist = [y_hist;y_next];
    fi_hist = [fi_hist; fi_next];
  
    fprintf("Iter: %d | Obj: %4.2f | Max_Cons: %1.2f \n", iter,y_next,max(fi_next));
    
    if norm(x_next-x_current)<=epsl
        disp('converged');
        break;
    end

end

% plot figure
plot_figure(x_hist,y_hist,fi_hist);
 

%% Auxiliary functions

% Define objective function
function y = obj_fun(x)
    y = (x(1)-2.7)^2+0.5*(x(2)-0.5)^2-5;
end

% Define constraint function
function fi = fi_fun(x)
    fi = 1.5*sin(x(1))-x(2)-0.1;
end

% Plot figures
function []=plot_figure(x_hist,y_hist,fi_hist)
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
    
    z = linspace(-3,3.5,100);
    f_1 = 1.5*sin(z)-0.1;
    plot(z,f_1,'r-',LineWidth=1);
    
    scatter(x_hist(1,1),x_hist(1,2),100,'ro','filled','MarkerFaceColor','r')
    scatter(2.7,0.5,100,'gp','filled','MarkerFaceColor','g')
    
    title('Optimization trajectory')
    xlabel('$x_1$','Interpreter','latex')
    ylabel('$x_2$','Interpreter','latex')
    legend('e-SLS','Constraint 1','Start','Optimum','Location','northwest')
    xlim([-3,3])
    ylim([-2,2])
end