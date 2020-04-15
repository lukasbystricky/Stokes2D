close all
clearvars
clc

% create input structure
input_params = default_input_params_periodic('curved_channel_demo');

% modify structure as needed, or add additional problem-dependent params
input_params.box_size = [5,1.5];
input_params.panels = 10;
input.pressure_drop_x = 1;

% set this up as a test problem. Shear flow boundary conditions are applied
% along with zero pressure gradient. Then the exact solution (shear flow)
% is known, so we can get an idea of how many panels are sufficient
input_params.test = 1;
input_params.plot_domain = 1;

input_params.n_periods_top = 2;
input_params.n_periods_bottom = 1;
input_params.amplitude_top = 0.1;
input_params.amplitude_bottom = 0.05;
input_params.h = 0.5;

problem = curved_pipe_periodic(input_params);

% solve the problem
solution = solve_stokes(problem);

% display solution

[Uc, Vc, X, Y,U] = evaluate_velocity(solution, 100);


h = figure();
if input_params.test
    
    plot_domain(problem, h);
    hold on
    
    exact_solution = @(x,y) y;
    contourf(X,Y, log10(abs(Uc - exact_solution(X,Y))+eps));
    colorbar
    axis equal
    title('log_{10}(error)');
    
else
   
    subplot(2,1,1);
    
    plot_domain(problem, h);
    hold on
    contourf(X,Y,Uc);
    colorbar
    axis equal
    title('U_1');
    
    subplot(2,1,2);
    
    plot_domain(problem, h);
    hold on    
    contourf(X,Y,Vc);
    colorbar
    axis equal
    title('U_2');
end



