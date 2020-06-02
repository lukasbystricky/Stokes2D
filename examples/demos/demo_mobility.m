% Solves the mobility problem for multiple particles in an unbounded
% domain. Compute the density q, as well a translational and rotational
% velocity. NOTE: THIS SIMULATION DOES NOT TIMESTEP. 

close all
clearvars
clc

n_high = 100;
n_low = 10;
h = 2*pi/n_low;
d = 0.05*h;

% create input structure
input_params = default_input_params('mobility', 0);

% modify structure as needed, or add additional problem-dependent params
input_params.plot_domain = 0;

% set up radii and centers, centers given as complex numbers (x+iy)
input_params.radii = [0.5, 0.5];
input_params.centers = [0, 1i*(1 + d)];

%% create high resolution solution
input_params.panels = n_high;

problem = particles(input_params);

% solve the problem
solution_high_res = solve_stokes(problem);

%% create low resolution solution
input_params.panels = n_low;

problem = particles(input_params);

% solve the problem
solution_low_res = solve_stokes(problem);

%% compute error in velocity and plot density
u_error = abs(solution_high_res.u_trans - solution_low_res.u_trans);

figure()
subplot(2,1,1)
plot(solution_high_res.problem.domain.theta(1:end/2), solution_high_res.q(1:end/2,:));
hold on
plot(solution_low_res.problem.domain.theta(1:end/2), solution_low_res.q(1:end/2,:));
title('particle 1');

subplot(2,1,2)
plot(solution_high_res.problem.domain.theta(end/2+1:end), solution_high_res.q(end/2+1:end,:));
hold on
plot(solution_low_res.problem.domain.theta(end/2+1:end), solution_low_res.q(end/2+1:end,:));
title('particle 2');