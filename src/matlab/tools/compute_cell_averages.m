function [u_avg, p_avg, p_grad_avg] = compute_cell_averages(solution, Nx, Ny)

u_avg = zeros(Nx*Ny,2);
p_avg = zeros(Nx*Ny,1);
p_grad_avg = zeros(Nx*Ny,2);

Lx = solution.problem.domain.Lx;
Ly = solution.problem.domain.Ly;
hx = Lx / Nx;
hy = Ly / Ny;

% compute velocity and velocity gradient on surface

% compute integral over boundary
[Ux, Uy, Vx, Vy] = evaluate_velocity_gradient_on_surface(solution);
[U, V] = evaluate_velocity_on_surface(solution);
[Px,Py] = evaluate_pressure_gradient_on_surface(solution);
P = evaluate_pressure_on_surface(solution);

% domain ranges over [-Lx/2, Lx/2]x[-Ly/2, Ly/2]
for i = 1:Nx
    for j = 1:Ny
        xmin = -Lx/2 + hx*(i-1);
        xmax = xmin + hx;
        ymin = -Ly/2 + hy*(j-1);
        ymax = ymin + hy;
        
        u_avg(Nx*(i-1) + j,:) = compute_average_velocity(solution, ...
                xmin, xmax, ymin, ymax, 400, 400, U, V, Ux, Uy, Vx, Vy) / (hx*hy);   
          
        p_avg(Nx*(i-1) + j) = compute_average_pressure(solution, ...
                xmin, xmax, ymin, ymax, 400, 400, P, Px, Py) / (hx*hy);  
            
        p_grad_avg(Nx*(i-1) + j,:) = compute_average_pressure_gradient(solution, ...
                xmin, xmax, ymin, ymax, 400, 400, P) / (hx*hy);        
    end
end