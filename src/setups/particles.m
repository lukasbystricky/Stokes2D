function problem = particles(input_params)
% Setup file a Couette apparatus, the outer wall is spinning, while all 
% other walls are fixed

problem.plot_domain = input_params.plot_domain;
problem.periodic = 0;
problem.name = input_params.name;
problem.gmres_tol = input_params.gmres_tol;
problem.stresslet_id_test = input_params.stresslet_id_test;
problem.resistance = false;

% discretize domain
radii = input_params.radii;
centers = input_params.centers;

particles = cell(length(radii));

% the first wall is always the bounding wall and has a normal pointing in
% the other direction to the other walls
for i = 1:length(radii)
    particles{i} =  @(T) geometry_circle(radii(i),centers(i),T);
end

if length(input_params.panels) > 1
    problem.panels = input_params.panels;
else
    problem.panels = input_params.panels*ones(length(particles),1);
end

problem.domain = discretize_domain(particles, problem.panels, centers);

% set up a separate domain for each particle, this will be used in the
% special quadrature
problem.domain.particles = cell(length(radii),1);
for i = 1:length(particles)
    problem.domain.particles{i} = discretize_domain({particles{i}}, problem.panels(i), centers(i));
end

% this is actually the background flow, which we will be extensional flow
problem.boundary_conditions = @(z) real(z) - 1i*imag(z);

if problem.plot_domain
    plot_domain(problem);
    disp('Simulation paused. Press any key to continue...');
    pause;
    disp('Continuing...');
end


