% Compute the permeability of an doubly-periodic array of circles. The 
% general procedure is outline for example in 
% https://arxiv.org/abs/1906.06884. For a reference cell containing a 
% single circle, with volume fraction phi=0.4, we can compare to the 
% results in that paper. There they get K = (5.671e-4, 0; 0, 5.671e-4). 

% NB: To do the volume meshing you will need the Matlab PDE toolbox

close all
clearvars
clc

% create input structure
input_params = default_input_params('permeability_demo', 1);

% modify structure as needed, or add additional problem-dependent params
input_params.box_size = [1,1];
input_params.panels = 10;
input_params.plot_domain = 0;

% stesslet identity has some issues if a boundary crosses over a periodic
% cell, so don't apply it here since quadrature points for volume integral
% should be inside the domain anyways
input_params.stresslet_id_test = @(t) zeros(size(t)); 

% volume fraction
phi = 0.4;

% set up radii and centers, centers given as complex numbers (x+iy)
input_params.radii = input_params.box_size(1)*sqrt((1-phi)/pi);

% for a single circle in a doubly periodic domain it doesn't matter where
% the center is, we've placed in at the corner of the reference cell to
% test the special quadrature, moving the center should not change the
% results
input_params.centers = 0;

K = zeros(2,2);

% solve for K_{11} and K_{12} by imposing pressure gradient in x direction
input.pressure_drop_x = 1;

problem = circles_periodic(input_params);

% solve the problem
solution = solve_stokes(problem);

figure(1);
% compute volume integral of u
[K(1,1), K(1,2)] = compute_volume_integral_circles(input_params.centers,...
    input_params.radii, solution, 'hmax', 0.02, 'N', 5);

% solve for K_{21} and K_{22} by imposing pressure gradient in y direction
input_params.pressure_drop_x = 0;
input_params.pressure_drop_y = 1;

problem = circles_periodic(input_params);

% solve the problem
solution = solve_stokes(problem);

figure(1);
% compute volume integral of u
[K(2,1), K(2,2)] = compute_volume_integral_circles(input_params.centers,...
    input_params.radii, solution, 'hmax', 0.02, 'N', 5);


fprintf("K_{11} = %3.3e\n", K(1,1));
fprintf("K_{12} = %3.3e\n", K(1,2));
fprintf("K_{21} = %3.3e\n", K(2,1));
fprintf("K_{22} = %3.3e\n", K(2,2));
