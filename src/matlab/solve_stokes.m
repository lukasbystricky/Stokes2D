function solution = solve_stokes(problem)

solution.problem = problem;

% set up right hand side
z = problem.domain.z;
nsrc = length(z);

rhs = [real(problem.boundary_conditions(z));... 
        imag(problem.boundary_conditions(z))];

if problem.periodic
    % add pressure constraint
    rhs = [rhs; problem.pressure_gradient_x; problem.pressure_gradient_y];

    X = gmres(@(x) matvec_combined(x, problem.domain, problem.eta), rhs, [], ...
                problem.gmres_tol, length(rhs));
else
    % add rows for net force and torque on inner walls
    nwalls = size(problem.domain.wall_indices, 1);    
    
    if problem.resistance 
        rhs = [-rhs; zeros(3*(nwalls-1), 1)];
        
        X = gmres(@(x) matvec_double_layer_resistance(x, problem.domain), rhs, [], ...
                problem.gmres_tol, length(rhs));
    else
        rhs = [rhs; zeros(3*(nwalls), 1)];
        
        X = gmres(@(x) matvec_double_layer_mobility(x, problem.domain), rhs, [], ...
                problem.gmres_tol, length(rhs));
    end
end
                        
solution.q = [X(1:nsrc), X(nsrc+1:2*nsrc)];

if problem.periodic
    % extract average velocity
    solution.u_avg = [X(end-1);X(end)];
else
   % extract net force and torques
   if problem.resistance
       n_inner_walls = length(problem.domain.centers) - 1;
       
       solution.forces = X(2*nsrc+1:2*nsrc+n_inner_walls) + ...
                1i*X(2*nsrc+n_inner_walls+1:2*nsrc+2*n_inner_walls);
       solution.torques = X(2*nsrc+2*n_inner_walls+1:end);   
   else % extract translational and rotational velocities
       n_particles = length(problem.domain.centers);
       
       solution.u_trans = X(2*nsrc+1:2*nsrc+n_particles) + ...
                1i*X(2*nsrc+n_particles+1:2*nsrc+2*n_particles);
       solution.omega = X(2*nsrc+2*n_particles+1:end);   
       
       solution.forces = zeros(n_particles,1);
       solution.torques = zeros(n_particles,1);
       
   end
end
    

