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
input_params.box_size = [5,1.5];
input_params.panels = 10;
 
problem = flat_pipe_periodic(input_params);

% solve the problem
solution = solve_stokes(problem);

% display solution
exact_solution = @(x,y) -problem.pressure_gradient_x * y.*(y-1)/2;

[Uc, ~, X, Y] = evaluate_velocity(solution, 100);

subplot(2,1,1)
contourf(X,Y,Uc);
colorbar
axis equal
title('U_1');

subplot(2,1,2)
contourf(X,Y, log10(abs(Uc - exact_solution(X,Y))+eps));
colorbar
axis equal
title('log_{10}(error)');



