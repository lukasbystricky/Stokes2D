%% Stress on-surface test
close all
clearvars
clc

% create input structure
input_params = default_input_params('periodic_starfish', 1);

% modify structure as needed
input_params.box_size = [5,3];
input_params.panels = 80;
input_params.centers = 0;
input_params.radii = 1;
input_params.plot_domain = 0;

problem = starfish_periodic(input_params);
problem.eta = 1;

% create solution structure
solution = solve_stokes(problem);
solution.local_indices = 1:length(solution.q);

% target points
theta = linspace(0,2*pi,1000) - 0.005*1i;
ztar = problem.walls{1}(theta);

%% evaluate stress both on and near surface
solution.problem.eta = 1;
[sxx_surf, sxy_surf, syx_surf, syy_surf] = evaluate_stress_on_surface(solution, solution, 'fluid');
[sxx_near, sxy_near, syx_near, syy_near] = evaluate_stress(solution, real(ztar), imag(ztar));

%% plot
figure;
subplot(2,2,1);
plot(problem.domain.theta, sxx_surf)
hold on
plot(real(theta), sxx_near);
legend('surface', 'near');
title('sxx');
xlim([0,2*pi]);
xlabel('theta');
ylabel('stress');

subplot(2,2,2);
plot(problem.domain.theta, sxy_surf)
hold on
plot(real(theta), sxy_near);
legend('surface', 'near');
title('sxy');
xlim([0,2*pi]);
xlabel('theta');
ylabel('stress');

subplot(2,2,3);
plot(problem.domain.theta, syx_surf)
hold on
plot(real(theta), syx_near);
legend('surface', 'near');
title('syx');
xlim([0,2*pi]);
xlabel('theta');
ylabel('stress');

subplot(2,2,4);
plot(problem.domain.theta, syy_surf)
hold on
plot(real(theta), syy_near);
legend('surface', 'near');
title('syy');
xlim([0,2*pi]);
xlabel('theta');
ylabel('stress');
