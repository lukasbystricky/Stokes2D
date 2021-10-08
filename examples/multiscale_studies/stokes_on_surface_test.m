% Compute pressure, pressure gradient, and velocity gradient on the surface
% when the boundary is approached from the fluid. The on-surface quantities
% are compared with them found using near surface evaluation.
close all
clearvars
clc

% create input structure
input_params = default_input_params('permeability_demo', 1);

% modify structure as needed, or add additional problem-dependent params
input_params.box_size = [5,3];
input_params.panels = 40;
input_params.plot_domain = 1;
input_params.eta = 1;

% prescribe pressure drop
input_params.pressure_drop_x = 1;
input_params.pressure_drop_y = 0;

problem = starfish_periodic(input_params);
% input_params.radii = 0.2;
% input_params.centers = 0;
% 
% problem = circles_periodic(input_params);

% create solution structure
solution = solve_stokes(problem);
solution.local_indices = 1:length(solution.q);

%% evaluate pressure both on and near surface
Psurface = evaluate_pressure_on_surface(solution, solution, 'fluid');
[Uxsurf, Uysurf, Vxsurf, Vysurf] = evaluate_velocity_gradient_on_surface(solution, solution, 'fluid');
[Pxsurf, Pysurf] = evaluate_pressure_gradient_on_surface(solution, solution, 'fluid');

[P,~, X, Y] = evaluate_pressure(solution, 50);
[Ux, Uy, Vx, Vy] = evaluate_velocity_gradient(solution, 50);
theta = linspace(0,2*pi,1000) - 0.005*1i;
% r = [input_params.radii + 0.001, input_params.radii - 0.001];
% x = [r(1)*cos(theta), r(2)*cos(theta)];
% y = [r(1)*sin(theta), r(2)*sin(theta)];
ztar = problem.walls{1}(theta);

solution.problem.stresslet_id_test = @(t) zeros(size(t));
Pnear = evaluate_pressure(solution, real(ztar), imag(ztar));
[Uxnear, Uynear, Vxnear, Vynear] = evaluate_velocity_gradient(solution, real(ztar), imag(ztar));
[dPxnear, dPynear] = evaluate_pressure_gradient(solution, real(ztar), imag(ztar));

%% compute average velocity
Lx = problem.domain.Lx;
Ly = problem.domain.Ly;
u_avg = compute_average_velocity(solution, solution, -Lx/2, Lx/2, -Ly/2, Ly/2, 20, 20);

%% plot
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
plot(real(theta), Pnear);

h=figure();
subplot(2,2,1);
contourf(X,Y,Ux);
axis equal; colorbar;
subplot(2,2,2);
contourf(X,Y,Uy);
axis equal; colorbar;
subplot(2,2,3);
contourf(X,Y,Vx);
axis equal; colorbar;
subplot(2,2,4);
contourf(X,Y,Vy);
axis equal; colorbar;

h=figure();
subplot(2,2,1)
plot(problem.domain.theta, Uxsurf);
hold on
plot(real(theta), Uxnear);
legend('Uxsurf','Uxnear');
title('Ux');

subplot(2,2,2)
plot(problem.domain.theta, Uysurf);
hold on
plot(real(theta), Uynear);

subplot(2,2,3)
plot(problem.domain.theta, Vxsurf);
hold on
plot(real(theta), Vxnear);

subplot(2,2,4)
plot(problem.domain.theta, Vysurf);
hold on
plot(real(theta), Vynear);

figure();
subplot(2,1,1);
plot(problem.domain.theta, Pxsurf);
hold on
plot(real(theta), dPxnear);

subplot(2,1,2);
plot(problem.domain.theta, Pysurf);
hold on
plot(real(theta), dPynear);