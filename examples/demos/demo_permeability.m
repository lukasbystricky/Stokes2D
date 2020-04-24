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
input_params.panels = 20;
input_params.plot_domain = 0;

% void fraction
phi = 0.4;

% set up radii and centers, centers given as complex numbers (x+iy)
input_params.radii = input_params.box_size(1)*sqrt((1-phi)/pi);

% for a single circle in a doubly periodic domain it doesn't matter where
% the center is
input_params.centers = 0.2 + 2.3*1i;

K = zeros(2,2);

%% solve for K_{11} and K_{12} by imposing pressure gradient in x direction
input.pressure_drop_x = 1;

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

fprintf("K_{11} = %3.3e\n", K(1,1));
fprintf("K_{12} = %3.3e\n", K(1,2));
fprintf("K_{21} = %3.3e\n", K(2,1));
fprintf("K_{22} = %3.3e\n", K(2,2));
