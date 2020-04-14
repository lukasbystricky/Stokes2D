function problem = periodic_pipe_setup(input_params)
% Setup file for a simple periodic pipe. The exact solution should be the 
% standard Pouseille flow.        

problem.Lx = input_params.box_size(1);
problem.Ly = input_params.box_size(2);

problem.plot_domain = input_params.plot_domain;
problem.periodic = 1;
problem.pressure_gradient_x = input_params.pressure_drop_x / problem.Lx;
problem.pressure_gradient_y = input_params.pressure_drop_y / problem.Ly;
problem.panels = input_params.panels;
problem.name = input_params.name;
problem.gmres_tol = input_params.gmres_tol;
problem.eta = input_params.eta;

% discretize domain
% top wall
walls{1} = @(T) geometry_periodic_channel(@(t) ones(size(t)), ...
            @(t) zeros(size(t)), @(t) zeros(size(t)),T, 1,problem.Lx);
     
% bottom wall            
walls{2} = @(T) geometry_periodic_channel(@(t) 0*ones(size(t)), ...
            @(t) zeros(size(t)), @(t) zeros(size(t)),T, -1, problem.Lx);
            
problem.domain = discretize_domain_periodic(walls, problem.panels,...
            problem.Lx, problem.Ly);

% specify boundary conditions, here we have simple no-slip
problem.boundary_conditions = @(z) 0*imag(z);

if problem.plot_domain
    plot_domain(problem);
    disp('Simulation paused. Press any key to continue...');
    pause;
end


