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
problem.alpha = input_params.alpha;

% discretize domain
% top wall
walls{1} = @(T) geometry_periodic_channel(@(t) problem.h*ones(size(t)), ...
            @(t) zeros(size(t)), @(t) zeros(size(t)),T, 1,problem.Lx);
     
% bottom wall            
walls{2} = @(T) geometry_periodic_channel(@(t) -problem.h*ones(size(t)), ...
            @(t) zeros(size(t)), @(t) zeros(size(t)),T, -1, problem.Lx);

% periodic sine as wall
% top wall
% A = 5*1e-2;
% walls{1} = @(T) geometry_periodic_channel(@(t) 0.5+A*sin(2*pi*t), ...
%             @(t) 2*pi*A*cos(2*pi*t), @(t) -4*pi*pi*A*sin(2*pi*t),T, 1,problem.Lx);
%      
% % bottom wall            
% walls{2} = @(T) geometry_periodic_channel(@(t) -0.5+A*sin(2*pi*t), ...
%             @(t) 2*pi*A*cos(2*pi*t), @(t) -4*pi*pi*A*sin(2*pi*t),T, -1, problem.Lx);
        
if length(input_params.panels) > 1
    problem.panels = input_params.panels;
else
    problem.panels = input_params.panels*ones(length(walls),1);
end

problem.domain = discretize_domain(walls, problem.panels,...
            problem.Lx, problem.Ly);

% specify boundary conditions, here we have simple no-slip
problem.boundary_conditions = @(z) 0*imag(z);

if problem.plot_domain
    plot_domain(problem);
    disp('Simulation paused. Press any key to continue...');
    pause;
end


