# Safe_Line_Search_Optimization

This is a model-free (black-box) safe optimization algorithm. It does not require the mathematical expressions of the objective and constraint functions of the optimization problem. Utilizing the lipschitz and smoothness assumptions about the unknown functions, this algorithm takes the inputs and outputs of the underlying functions to determine the next iteration point. It is guaranteed that no evaluation during optimization shall be made. The algorithm has two versions: 
* Algorithm for exact measurements (e-SLS);
* Algorithm for noisy measurements (n-SLS);

# Algorithm outline
The algorithm consists of the following steps:
1. Estimate the gradients of the unknown functions;
2. Determine a search direction;
3. Compute a safe set for the next iterations;
4. Project the search direction if near to a boundary;
5. Compute a safe step length;
6. Set the next iteration point;
7. Evaluate functions;
8. Check convergence condition.

# Usage
The esls.m function in the interface directory is a callable function of the algorithm. It takes the following inputs:
the current iteration point, objective function handle, constraint function handle, algorithm parameters, and return:
the next iteration point, lagrangian vector, convergence index. A usage example of the function is the esls_example.m

# Examples
The files e_SLS_case1.m and n_SLS_case1.m ar 2D examples of the algorithm for the next optimizaitoin problem with exact measurements and noisy measurements respectively:  
&emsp; min (x_1-2.7)^2+0.5(x_2-0.5)-5  
&emsp; s.t. &ensp;  x1<=2.7  
&emsp; &emsp;&emsp;  x2>=-5
