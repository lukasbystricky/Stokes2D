function problem = boundary_layer_test(input_params)
% Setup file for a doubly periodic array of circles       

problem.Lx = input_params.box_size(1);
problem.Ly = input_params.box_size(2);

problem.plot_domain = input_params.plot_domain;
problem.periodic = 1;
problem.pressure_gradient_x = input_params.pressure_drop_x / problem.Lx;
problem.pressure_gradient_y = input_params.pressure_drop_y / problem.Ly;
problem.name = input_params.name;
problem.gmres_tol = input_params.gmres_tol;
problem.eta = input_params.eta;
problem.stresslet_id_test = input_params.stresslet_id_test;

% discretize domain
radii = input_params.radii;
centers = input_params.centers;

walls = cell(length(radii) + 2,1);

for i = 1:length(radii)
   walls{i} =  @(T) geometry_circle(radii(i),centers(i),T);
end

% top wall
walls{end-1} = @(T) geometry_periodic_channel(@(t) (problem.Ly/2 - 1)*ones(size(t)), ...
            @(t) zeros(size(t)), @(t) zeros(size(t)),T, 1,problem.Lx);
     
% bottom wall            
walls{end} = @(T) geometry_periodic_channel(@(t) -problem.Ly/2*ones(size(t)), ...
            @(t) zeros(size(t)), @(t) zeros(size(t)),T, -1, problem.Lx);

if length(input_params.panels) > 1
    problem.panels = input_params.panels;
else
    problem.panels = input_params.panels*ones(length(walls),1);
end

problem.domain = discretize_domain(walls, problem.panels,...
            problem.Lx, problem.Ly);

% specify boundary conditions
if input_params.test
    % test, no pressure gradient, apply shear flow
    problem.pressure_gradient_x = 0;
    problem.pressure_gradient_y = 0;
    problem.boundary_conditions = @(z) imag(z);
else
    % no-slip
    problem.boundary_conditions = @(z) 0*imag(z);
    
    % shear on top boundary
%     problem.boundary_conditions = @(z) 0.1*(abs(imag(z) - (problem.Ly/2-1)) < 1e-2);
%     problem.pressure_gradient_x = 0;
end

if problem.plot_domain
    plot_domain(problem);
    disp('Simulation paused. Press any key to continue...');
    pause;
    disp('Continuing...');
end


