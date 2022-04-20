function [x_next,lambda,converged,G0] = esls(x_current,x_last,gt_last,obj_hand,fi_hand,options)
%This function runs the safe line-search optimization algorithm.
%
%   [x_next,lambda,converged] = esls(x_current,obj_hand,fi_hand,options) 
%   
%   Parameters
%   ------------------
%   x_current: array
%       Decision point at the current iteration
%   x_last: array or nan
%       Decision point at the last iteration
%   gt_last: array or nan
%       Gradient of the objective function at the last iteration
%   obj_hand: function handle
%       Objective function_handle, example: y_current=obj_fun(x_current)
%   fi_hand: function handle
%       Constraint function_handle, example: fi_current=fi_fun(x_current)
%   options: struct
%       Algorithm parameters:
%           L:    real
%                 Lipschitz constant
%           M:    real
%           mu:   real
%                 Gradient estimation deviation upperbound
%           h:    real
%                 Safety threshold
%           epsl: real
%                 Convergence condition
%           rho:  real;
%                 Update rate of step length selection
%           c:    real
%                 Small constant in step length selection
%       example: options.L = 1, options.M = 1
%   ------------------
%   Return:
%       x_next: array
%           Decision point at the next iteration
%       lambda: array
%           Lagrangian multiplier vector
%       converged: int
%           convergence index, 1-converged, 0-not converged


    y_current = obj_hand(x_current);
    fi_current = fi_hand(x_current);
    
    d = length(x_current);
    m = length(fi_current);
    L = options.L;
    M = options.M;
    mu = options.mu;
    h = options.h;
    epsl = options.epsl;
    rho = options.rho;
    c = options.c;
    use_newton_direction = options.newton_direction;
    
    lambda = zeros(1,m);
    H = eye(d);

    vk = min(2*mu/(sqrt(d)*M), min(-fi_current)/(2*L));   % Compute step length for finite difference
    
    X_current = repmat(x_current,d,1);
    VT = eye(d)*vk;
    X_current = X_current+VT;
    Y_current = [];
    FI_current = [];
    for i = 1:d
        y_e = obj_hand(X_current(i,:));
        fi_e = fi_hand(X_current(i,:));
        Y_current = [Y_current; y_e];
        FI_current = [FI_current; fi_e];
    end
    
    G0 = (Y_current - y_current)/vk;        % Compute gradient estimator for objective
    GI = (FI_current - fi_current)/vk;      % Compute gradient estimator for constraints
    
    % Define the search direction
    if use_newton_direction
        if prod(~isnan(x_last))&&prod(~isnan(gt_last))
            x_diff = x_current'-x_last';
            gt_diff = G0 - gt_last';
            psi = 1/(gt_diff'*x_diff);
            H = (eye(d) - psi*x_diff*gt_diff')*H*(eye(d)-psi*gt_diff*x_diff') + psi*(x_diff*x_diff');
        end

        p = -H*G0;
    else
        p = -G0;                            
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
    
    y_k1 = obj_hand(x_current+alpha*p');
    fi_k1 = fi_hand(x_current+alpha*p');
    
    % Select safe step length
    while y_k1 > y_k + c*alpha*(G0'*p)
        alpha = alpha* rho;
    
        y_k1 = obj_hand(x_current+alpha*p');
        fi_k1 = fi_hand(x_current+alpha*p');
    end
    
    % Update iteration
    x_next = x_current + alpha*p';
    
    if norm(x_next-x_current)<=epsl
        converged = 1;
    else
        converged = 0;
    end

end