function solution = solve_stokes(problem)

solution.domain = problem.domain;

% set up right hand side
z = problem.domain.z;

rhs = [real(problem.boundary_conditions(z));... 
        imag(problem.boundary_conditions(z))];

% add pressure constraint
rhs = [rhs; problem.pressure_gradient_x; problem.pressure_gradient_y];

X = gmres(@(x) matvec_combined(x, problem.domain, problem.eta), rhs, [], ...
                problem.gmres_tol, length(rhs));
                        
solution.q = [X(1:length(z)), X(length(z)+1:2*length(z))];
solution.u_avg = [X(end-1);X(end)];
solution.eta = problem.eta;            
    

