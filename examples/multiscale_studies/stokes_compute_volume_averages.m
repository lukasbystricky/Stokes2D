% Compute the permeability of an doubly-periodic array of obstacles. The 
% general procedure is outline for example in 
% https://arxiv.org/abs/1906.06884. For a reference cell containing a 
% single circle, with volume fraction phi=0.4, we can compare to the 
% results in that paper. There they get K = (5.671e-4, 0; 0, 5.671e-4). 

close all
%clearvars
clc

% create input structure
input_params = default_input_params('permeability_demo', 1);

% modify structure as needed, or add additional problem-dependent params
input_params.box_size = [1,1];
input_params.panels = 5;
input_params.plot_domain = 0;
input_params.eta = 1;

% void fraction
phi = 0.8;

% set up radii and centers, centers given as complex numbers (x+iy)
input_params.radii = input_params.box_size(1)*sqrt((1-phi)/pi);

% for a single circle in a doubly periodic domain it doesn't matter where
% the center is
input_params.centers = 0;

K = zeros(2,2);

%% solve for K_{11} and K_{12} by imposing pressure gradient in x direction
input_params.pressure_drop_x = 1;
input_params.pressure_drop_y = 0;

problem = circles_periodic(input_params);

% solve the problem
solution = solve_stokes(problem);

% average velocity is computed already!
K(:,1) = solution.u_avg;

%% solve for K_{21} and K_{22} by imposing pressure gradient in y direction
input_params.pressure_drop_x = 0;
input_params.pressure_drop_y = 1;

problem = circles_periodic(input_params);

% solve the problem
solution = solve_stokes(problem);

K(:,2) = solution.u_avg;

fprintf("SINGLE CIRCLE:\n");
fprintf("K_{11} = %3.3e\n", K(1,1));
fprintf("K_{12} = %3.3e\n", K(1,2));
fprintf("K_{21} = %3.3e\n", K(2,1));
fprintf("K_{22} = %3.3e\n", K(2,2));

Psurface = evaluate_pressure_on_surface(solution);
[P,~, X, Y] = evaluate_pressure(solution, 100);

theta = linspace(0,2*pi);
r = [input_params.radii + 0.001, input_params.radii - 0.001];
x = [r(1)*cos(theta), r(2)*cos(theta)];
y = [r(1)*sin(theta), r(2)*sin(theta)];

solution.problem.stresslet_id_test = @(t) zeros(size(t));
Pnear = evaluate_pressure(solution, x, y);

h=figure();
subplot(1,2,1);
contourf(X,Y,P);
hold on
plot_domain(problem, h);
axis equal;
colorbar;
 
subplot(1,2,2);
plot(problem.domain.theta, Psurface)
hold on
plot(theta, Pnear(1:end/2));

