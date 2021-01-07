function [u_avg, p_avg, p_grad_avg, u_grad_avg] = compute_cell_averages(solution, x, y, hx, hy, bodies)

Nx = length(x);
Ny = length(y);
Ncells = Nx*Ny;
u_avg = zeros(Ncells,2);
p_avg = zeros(Ncells);
p_grad_avg = zeros(Ncells,2);
u_grad_avg = zeros(Ncells,2, 2);

% compute velocity and velocity gradient on surface

% compute integral over boundary
% [Ux, Uy, Vx, Vy] = evaluate_velocity_gradient_on_surface(solution, 1);
% [U, V] = evaluate_velocity_on_surface(solution,1);
% [Px,Py] = evaluate_pressure_gradient_on_surface(solution,1);
% P = evaluate_pressure_on_surface(solution,1);

% domain ranges over [-Lx/2, Lx/2]x[-Ly/2, Ly/2]
for i = 1:length(x)
    for j = 1:length(y)        
            xmin = x(i) - hx/2;
            xmax = xmin + hx;
            ymin = y(j);
            ymax = ymin + hy;
            
            u_avg(Nx*(i-1) + j,:) = compute_average_velocity(solution, ...
                xmin, xmax, ymin, ymax, 400, 400, bodies) / (hx*hy);
            
            p_avg(Nx*(i-1) + j) = compute_average_pressure(solution, ...
                xmin, xmax, ymin, ymax, 400, 400, bodies) / (hx*hy);
            
            p_grad_avg(Nx*(i-1) + j,:) = compute_average_pressure_gradient(solution, ...
                xmin, xmax, ymin, ymax, 400, 400) / (hx*hy);
            
            u_grad_avg(Nx*(i-1) + j,:, :) = compute_average_velocity_gradient(solution, ...
                xmin, xmax, ymin, ymax, 400, 400) / (hx*hy);
    end
end