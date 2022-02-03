function solution = solve_stokes(problem, rhs, varargin)
% default parameters
fmm = 1;
verbose = 0;

% read in optional input parameters
if nargin > 2
    % Go through all other input arguments and assign parameters
    jv = 1;
    while jv <= length(varargin)-1
       switch varargin{jv}
           case 'fmm'
               fmm = varargin{jv+1}; 
           case 'verbose'
               verbose = varargin{jv+1};
       end
       jv = jv + 2;
    end
end

solution.problem = problem;

% set up right hand side
z = problem.domain.z;
nsrc = length(z);

% rhs = [real(problem.boundary_conditions(z));... 
%         imag(problem.boundary_conditions(z))];

if problem.periodic
    % add pressure constraint
    %rhs = [rhs; problem.pressure_gradient_x; problem.pressure_gradient_y];
    
    if problem.slip
        X = gmres(@(x) matvec_combined_slip(x, problem), rhs, [], ...
            problem.gmres_tol, min(length(rhs),1000));
    else
        X = gmres(@(x) matvec_combined(x, problem.domain, problem.eta), rhs, [], ...
        problem.gmres_tol, min(length(rhs),500));
    end
else
    % add rows for net force and torque on inner walls
    nwalls = size(problem.domain.wall_indices, 1);    
    
    if problem.resistance
        rhs = [rhs; zeros(3*(nwalls-1), 1)];
        
        X = gmres(@(x) matvec_double_layer_resistance(x, problem.domain, fmm), rhs, [], ...
                problem.gmres_tol, length(rhs));
    else
        rhs = [-rhs; real(problem.forces);imag(problem.forces);problem.torques];
        
        [uS, uR] = completion_contribution(problem.domain.centers, z, ...
                    problem.forces, problem.torques);
                
        rhs(1:length(z)) = rhs(1:length(z)) + real(uS + uR);
        rhs(length(z)+1:2*length(z)) = rhs(length(z)+1:2*length(z)) + imag(uS + uR);
         
        X = gmres(@(x) matvec_double_layer_mobility(x, problem.domain, fmm), rhs, [], ...
                problem.gmres_tol, length(rhs));
    end
end

solution.q = [X(1:nsrc), X(nsrc+1:2*nsrc)];

if problem.periodic
    % extract average velocity
    solution.u_avg = [X(end-1);X(end)];
    
    % trim post-processed results by removing all values outside the domain
    solution.trim = 1;
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
       
       solution.forces = problem.forces;
       solution.torques = problem.torques;       
   end
end
    

