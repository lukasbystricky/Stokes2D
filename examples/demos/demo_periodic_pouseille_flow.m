% Basic demo script. Computes the velocity field for a pressure driven flow
% through a flat pipe, i.e. Pouseille flow. The exact solution is known and
% given by -y*(y-1)*dP/Lx/2, where dP/Lx is the imposed pressure gradient.
% Note that this problem is only periodic in the x direction, however we
% can still use the 2P periodic spectral Ewald code be embedding the 
% channel in a periodic box. We use a combined-layer formulation where the 
% pressure gradient is balanced by the single-layer potential. 

close all
clearvars
clc

% create input structure
input_params = default_input_params('pouseuille_demo', 1);

% modify structure as needed
input_params.box_size = [5,3];
input_params.panels = 5;
input_params.eta = 1;

problem = flat_pipe_periodic(input_params);

% solve the problem
solution = solve_stokes(problem);

% display solution
exact_solution_u = @(x,y) -problem.pressure_gradient_x * y.*(y-1)/2;
exact_solution_uy = @(x,y) -problem.pressure_gradient_x * (y - 0.5);
exact_solution_p = @(x,y) problem.pressure_gradient_x * x;
exact_solution_dP = @(x,y) problem.pressure_gradient_x;

x = linspace(0,5,20);
y = linspace(0,1,20);

[U, V, X, Y] = evaluate_velocity(solution, 100);
P = evaluate_pressure(solution, X, Y);
[Ux, Uy, Vx, Vy] = evaluate_velocity_gradient(solution, X, Y);
[Px, Py] = evaluate_pressure_gradient(solution, X, Y);

subplot(4,1,1)

contourf(X,Y, log10(abs((U + 1i*V - exact_solution_u(X,Y))./...
    max(max(abs(exact_solution_u(X,Y)))))+eps));
colorbar
caxis([-16,-8]);
axis equal
title('velocity: log_{10}(relative error)');

subplot(4,1,2)
contourf(X,Y,log10(abs(P - exact_solution_p(X,Y))+eps));
colorbar
caxis([-16,-8]);
axis equal
title('P: log_{10}(relative error)');

subplot(4,1,3)
contourf(X,Y,log10(abs(Uy - exact_solution_uy(X,Y))+eps));
colorbar
caxis([-16,-8]);
axis equal
title('u_y: log_{10}(relative error)');

subplot(4,1,4)
contourf(X,Y,log10(abs(Px + 1i*Py  - exact_solution_dP(X,Y))+eps));
colorbar
caxis([-16,-8]);
axis equal
title('nabla p: log_{10}(relative error)');


