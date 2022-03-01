function problem = couette(input_params)
% Setup file a Couette apparatus, the outer wall is spinning, while all 
% other walls are fixed       

problem.plot_domain = input_params.plot_domain;
problem.periodic = 0;
problem.name = input_params.name;
problem.gmres_tol = input_params.gmres_tol;
problem.stresslet_id_test = input_params.stresslet_id_test;
problem.resistance = true;
problem.slip = input_params.slip;
if problem.slip
    problem.alpha = input_params.alpha;
end

% discretize domain
radii = input_params.radii;
centers = input_params.centers;
omega = input_params.omega;

walls = cell(length(radii));

% the first wall is always the bounding wall and has a normal pointing in
% the other direction to the other walls
walls{1} =  @(T) geometry_circle(radii(1),centers(1),T, -1);

if length(radii) > 1
    for i = 2:length(radii)
        walls{i} =  @(T) geometry_circle(radii(i),centers(i),T);
    end
end

if length(input_params.panels) > 1
    problem.panels = input_params.panels;
else
    problem.panels = input_params.panels*ones(length(walls),1);
end

problem.domain = discretize_domain(walls, problem.panels, centers);

% zero on all walls except the outer, which is given by omega*[-y, x]
problem.boundary_conditions = ...
  @(z) omega*(abs(abs(z)-radii(2))<1e-12).*(-imag(z) + 1i*real(z));
% problem.boundary_conditions = ...
%   @(z) (-imag(z) + 1i*real(z))./(abs(z).^2);

if problem.plot_domain
    plot_domain(problem);
    disp('Simulation paused. Press any key to continue...');
    pause;
    disp('Continuing...');
end


