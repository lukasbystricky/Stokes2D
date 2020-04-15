function problem = curved_pipe_periodic(input_params)
% Setup file for a curved channel       

problem.Lx = input_params.box_size(1);
problem.Ly = input_params.box_size(2);

problem.plot_domain = input_params.plot_domain;
problem.periodic = 1;
problem.pressure_gradient_x = input_params.pressure_drop_x / problem.Lx;
problem.pressure_gradient_y = input_params.pressure_drop_y / problem.Ly;
problem.name = input_params.name;
problem.gmres_tol = input_params.gmres_tol;
problem.eta = input_params.eta;

k1 = 2*pi*input_params.n_periods_top/problem.Lx;
k2 = 2*pi*input_params.n_periods_bottom/problem.Lx;
a1 = input_params.amplitude_top;
a2 = input_params.amplitude_top;
h = input_params.h;

% discretize domain
% top wall
walls{1} = @(T) geometry_periodic_channel(@(t) a1*sin(k1*t) + h, ...
            @(t) a1*k1*cos(k1*t), @(t) -a1*k1^2*sin(k1*t),T, 1,problem.Lx);
     
% bottom wall            
walls{2} = @(T) geometry_periodic_channel(@(t) a2*sin(k2*t), ...
            @(t) a2*k2*cos(k2*t), @(t) -a2*k2^2*sin(k2*t),T, -1, problem.Lx);

if length(input_params.panels) > 1
    problem.panels = input_params.panels;
else
    problem.panels = input_params.panels*ones(length(walls),1);
end

problem.domain = discretize_domain_periodic(walls, problem.panels,...
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
end

if problem.plot_domain
    plot_domain(problem);
    disp('Simulation paused. Press any key to continue...');
    pause;
    disp('Continuing...');
end


