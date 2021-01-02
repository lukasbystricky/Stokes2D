function [u_avg, p_avg, p_grad_avg, u_grad_avg, x_cen, y_cen] = compute_cell_averages(solution, Nx, Ny, offsets)

n_offsets = length(offsets);

u_avg = zeros(Nx*Ny,2, n_offsets);
p_avg = zeros(Nx*Ny, n_offsets);
p_grad_avg = zeros(Nx*Ny,2, n_offsets);
u_grad_avg = zeros(Nx*Ny,2, 2, n_offsets);
x_cen = zeros(Nx*Ny, n_offsets);
y_cen = zeros(Nx*Ny, n_offsets);

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
for k = 1:n_offsets
    for i = 1:Nx
        for j = 1:Ny
            xmin = -Lx/2 + hx*(i-1);
            xmax = xmin + hx;
            ymin = -Ly/2 + hy*(j-1) + offsets(k);
            ymax = ymin + hy;
            
            x_cen(Nx*(i-1)+j, k) = (xmax + xmin)/2;
            y_cen(Nx*(i-1)+j, k) = (ymax + ymin)/2;
            
            u_avg(Nx*(i-1) + j,:, k) = compute_average_velocity(solution, ...
                xmin, xmax, ymin, ymax, 400, 400, U, V, Ux, Uy, Vx, Vy) / (hx*hy);
            
            p_avg(Nx*(i-1) + j, k) = compute_average_pressure(solution, ...
                xmin, xmax, ymin, ymax, 400, 400, P, Px, Py) / (hx*hy);
            
            p_grad_avg(Nx*(i-1) + j,:, k) = compute_average_pressure_gradient(solution, ...
                xmin, xmax, ymin, ymax, 400, 400, P) / (hx*hy);
            
            u_grad_avg(Nx*(i-1) + j,:, :, k) = compute_average_velocity_gradient(solution, ...
                xmin, xmax, ymin, ymax, 400, 400, P) / (hx*hy);
        end
    end
end