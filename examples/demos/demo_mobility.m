% Solves the mobility problem for multiple particles in an unbounded
% domain. Compute the density q, as well a translational and rotational
% velocity. NOTE: THIS SIMULATION DOES NOT TIMESTEP. 

close all
clearvars
clc

n_high = 200;
n_low = [10,20,40];
h = pi./n_low;

d = [10, 5, 1, 0.1, 0.01, 0.001]; 
u_error= zeros(length(n_low), length(d));

% create input structure
input_params = default_input_params('mobility', 0);

% modify structure as needed, or add additional problem-dependent params
input_params.plot_domain = 0;

% set up radii and centers, centers given as complex numbers (x+iy)
input_params.radii = [0.5, 0.5];

% set force and torque on each particle
input_params.forces = [1+1i; 0];
input_params.torques = [0; 2];

for i = 1:length(n_low)
    for j = 1:length(d)
        
        fprintf("n_pan = %d, d = %3.3g\n\n", n_low(i), d(j));
        
        input_params.centers = [0, 1i*(1 + h(i)*d(j))];
        
        %% create high resolution solution
        input_params.panels = n_high;
        
        problem = particles(input_params);
       
        
        % solve the problem
        solution_high_res = solve_stokes(problem);
        
        %% create low resolution solution
        input_params.panels = n_low(i);
        
        problem = particles(input_params);
        problem.forces = [1+1i; 0];
        problem.torques = [0; 2];
        
        % solve the problem
        solution_low_res = solve_stokes(problem);
        
        %% compute error in velocity and plot density
        u_error(i,j) = norm(abs(solution_high_res.u_trans - solution_low_res.u_trans));
        
        
    end
end

figure()
loglog(d, u_error)
xlabel("distance/h")
ylabel("error")

figure()
subplot(2,1,1)
plot(solution_high_res.problem.domain.theta(1:end/2), solution_high_res.q(1:end/2,:), '-o');
hold on
plot(solution_low_res.problem.domain.theta(1:end/2), solution_low_res.q(1:end/2,:), '-o');
title('particle 1');

subplot(2,1,2)
plot(solution_high_res.problem.domain.theta(end/2+1:end), solution_high_res.q(end/2+1:end,:), '-o');
hold on
plot(solution_low_res.problem.domain.theta(end/2+1:end), solution_low_res.q(end/2+1:end,:), '-o');
title('particle 2');