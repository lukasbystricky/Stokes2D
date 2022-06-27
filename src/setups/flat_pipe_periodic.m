function problem = flat_pipe_periodic(input_params)
% Setup file for a simple periodic pipe. The exact solution should be the 
% standard Pouseille flow.        

problem.Lx = input_params.box_size(1);
problem.Ly = input_params.box_size(2);

problem.plot_domain = input_params.plot_domain;
problem.periodic = 1;
problem.h = input_params.h;
problem.pressure_gradient_x = input_params.pressure_drop_x / problem.Lx;
problem.pressure_gradient_y = input_params.pressure_drop_y / problem.Ly;
problem.name = input_params.name;
problem.gmres_tol = input_params.gmres_tol;
problem.eta = input_params.eta;
problem.stresslet_id_test = input_params.stresslet_id_test;
problem.slip = input_params.slip;

if problem.slip
    problem.alpha = input_params.alpha;
else
    problem.alpha = 0;
end

% discretize domain
% top wall
walls{1} = @(T) geometry_periodic_channel(@(t) problem.h*ones(size(t)), ...
            @(t) zeros(size(t)), @(t) zeros(size(t)),T, 1,problem.Lx);
     
% bottom wall            
walls{2} = @(T) geometry_periodic_channel(@(t) -problem.h*ones(size(t)), ...
            @(t) zeros(size(t)), @(t) zeros(size(t)),T, -1, problem.Lx);
        
if length(input_params.panels) > 1
    problem.panels = input_params.panels;
else
    problem.panels = input_params.panels*ones(length(walls),1);
end

problem.domain = discretize_domain(walls, problem.panels,...
            problem.Lx, problem.Ly);

if isfield(input_params,'nbr_neighbor_pts')
    problem.domain.nbr_neighbor_pts = input_params.nbr_neighbor_pts;
else
    problem.domain.nbr_neighbor_pts = 4;
end

% specify boundary conditions, here we have simple no-slip
problem.boundary_conditions = @(z) 0*imag(z);

if problem.plot_domain
    plot_domain(problem);
    disp('Simulation paused. Press any key to continue...');
    pause;
end


