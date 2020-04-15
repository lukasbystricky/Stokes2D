close all
clearvars
clc

% create input structure
input_params = default_input_params_periodic('pouseuille_demo');

% modify structure as needed
input_params.box_size = [5,1.5];
input_params.panels = 10;
 
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



