function problem = curved_rectangle(input_params)
problem.plot_domain = input_params.plot_domain;
problem.periodic = 0;
problem.name = input_params.name;
problem.stresslet_id_test = input_params.stresslet_id_test;
problem.gmres_tol = input_params.gmres_tol;
problem.resistance = true;
problem.nsub = input_params.nsub;

% setups walls
p1 = input_params.corners(1);
p2 = input_params.corners(2);
p3 = input_params.corners(3);
p4 = input_params.corners(4);

walls = cell(4,1);

% top wall
walls{1} = @(T) geometry_periodic_channel(@(t) imag(p1)*ones(size(t)), ...
    @(t) zeros(size(t)), @(t) zeros(size(t)),T, 1, abs(p2-p1));
     
% right wall
walls{2} = @(T) vertical_line(T,real(p2),imag(p3),imag(p2),-1);

% bottom wall
A = 5*1e-2;
walls{3} = @(T) geometry_periodic_channel(@(t) imag(p4)+A*sin(2*pi*t), ...
    @(t) 2*pi*A*cos(2*pi*t), @(t) -4*pi*pi*A*sin(2*pi*t),T, -1, abs(p2-p1));

% walls{3} = @(T) geometry_periodic_channel(@(t) imag(p3)*ones(size(t)), ...
%     @(t) zeros(size(t)), @(t) zeros(size(t)),T, -1, abs(p2-p1));

% left wall
walls{4} = @(T) vertical_line(T,real(p4),imag(p4),imag(p1),1);

% walls{1} =  @(T) geometry_circle(1,0,T,-1);

% discretize domain
if length(input_params.panels) > 1
    problem.panels = input_params.panels;
else
    problem.panels = input_params.panels*ones(length(walls),1);
end

if problem.nsub > 0
    if any(problem.panels < 2)
        error('each boundary must consist of atleast two panels in order to dyadically subdivide');
    else
        problem.domain = discretize_domain_dyad_refine(walls, problem.panels, problem.nsub);
    end
else
    problem.domain = discretize_domain(walls, problem.panels);
end

% additional domain properties
problem.domain.n_inner_walls = 0;
problem.domain.inner_wall_indices = [];
problem.domain.outer_wall_indices = 1:length(problem.domain.z);

% set boundary conditions
problem.boundary_conditions = @(z) input_params.boundary_conditions(z);

% problem.boundary_conditions = @(z) 1*(abs(abs(z)-1)<1e-12).*(-imag(z) + 1i*real(z));

if problem.plot_domain
    plot_domain(problem);
    disp('Simulation paused. Press any key to continue...');
    pause;
    disp('Continuing...');
end
