% Periodic pressure driven flow through a channel with non-flat walls. Here
% we don't know the exact solution, but we can get an idea of the error by 
% imposing shear flow boundary conditions, and zero pressure gradient.

close all
clearvars
clc

% create input structure
input_params = default_input_params('curved_channel_demo', 1);

% modify structure as needed, or add additional problem-dependent params
input_params.box_size = [2,1.5];
input_params.panels = 50;
input_params.pressure_drop_x = 1;
input_params.eta = 1;

% set this up as a test problem. Shear flow boundary conditions are applied
% along with zero pressure gradient. Then the exact solution (shear flow)
% is known, so we can get an idea of how many panels are sufficient
input_params.test = 1;
input_params.plot_domain = 0;

% geometry information, periods, amplitude and some measure of spacing
% between walls
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
Pc  = evaluate_pressure(solution, X, Y);
[Uxc, Uyc, Vxc, Vyc, Ux, Uy, Vx, Vy] = evaluate_velocity_gradient(solution, X, Y);
[Px, Py]  = evaluate_pressure_gradient(solution, X, Y);

h = figure();
if input_params.test
    
    subplot(2,1,1)
    plot_domain(problem, h);
    hold on
    
    exact_solution = @(x,y) y;
    contourf(X,Y, log10(abs((Uc - exact_solution(X,Y))./...
                max(max(abs(exact_solution(X,Y)))))+eps));
    colorbar
    caxis([-16,-1]);
    axis equal
    title('velocity: log_{10}(relative error)');
    
    subplot(2,1,2)
    plot_domain(problem, h);
    hold on
    
    %exact_solution = @(x,y) mean(mean(Pc(~isnan(Pc)))); % constant pressure
    exact_solution = @(x,y) zeros(size(x)); % constant pressure
    contourf(X,Y, log10(abs((Px + 1i*Py - exact_solution(X,Y)))+eps));
   
    colorbar
    caxis([-16,-1]);
    axis equal
    title('pressure: log_{10}(relative error)');
    
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
    contourf(X,Y,Pc);
    colorbar
    axis equal
    title('Pressure');
end



