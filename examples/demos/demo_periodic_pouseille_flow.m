close all
clearvars
clc

% input strucutre:
% name: simulation name, data for the simulation will be saved somewhere
% panels: number of panels on each wall, vector of size 1 x number of walls
% plot_domain: flag specifying whether to plot the domain
% pressure_drop: applied pressure drop, specify only if domain is periodic
% box_size: size of the periodic box, a 1x2 vector [Lx, Ly]
% gmres_tol: GMRES tolerance
% eta: scaling parameter between single- and double-layers
input_params = struct('name', 'periodic_pipe',...
                      'panels', [10, 10],...
                      'plot_domain', 0,...
                      'pressure_drop_x', 1,...
                      'pressure_drop_y', 0,...
                      'box_size', [5, 1.5],...
                      'gmres_tol', 1e-12,...
                      'eta', 1);
 
problem = periodic_pipe_setup(input_params);

% solve the problem
solution = solve_stokes(problem);

% display solution
exact_solution = @(x,y) -problem.pressure_gradient_x * y.*(y-1)/2;

[Uc, ~, X, Y,U] = evaluate_velocity(solution, 100);

subplot(2,1,1)
contourf(X,Y,Uc);
colorbar
axis equal
title('U_1');

subplot(2,1,2)
contourf(X,Y, log10(abs(Uc - exact_solution(X,Y))+eps));
colorbar
axis equal
title('error');



