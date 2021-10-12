function [u_avg, p_avg, p_grad_avg, u_grad_avg] = compute_cell_averages(solution, x, y, hx, hy)

Nx = length(x);
Ny = length(y);
Ncells = Nx*Ny;
u_avg = zeros(Ncells,2);
p_avg = zeros(Ncells);
p_grad_avg = zeros(Ncells,2);
u_grad_avg = zeros(Ncells,2, 2);
domain = solution.problem.domain;
z = domain.z;

% cell ranges over [-Lx/2, Lx/2]x[-Ly/2, Ly/2]
for i = 1:length(x)
    for j = 1:length(y)        
            xmin = x(i) - hx/2;
            xmax = xmin + hx;
            ymin = y(j) - hy/2;
            ymax = ymin + hy;
            
            indices = [];
            wall_indices = [];
            bodies = [];
            wall_start = 1;
            n_bodies = size(domain.wall_indices,1);
            
            for k = 1:n_bodies
                indices_tmp = domain.wall_indices(k,1):domain.wall_indices(k,2);
                if min(imag(z(indices_tmp))) > ymin && max(imag(z(indices_tmp))) < ymax && ...
                        min(real(z(indices_tmp))) > xmin && max(real(z(indices_tmp))) < xmax
                    
                    indices = [indices, indices_tmp];
                    wall_indices = [wall_indices; wall_start, wall_start + length(indices_tmp)-1];
                    bodies = [bodies,k];
                    wall_start = wall_start + length(indices_tmp);
                end
            end
            
            solution_local = create_local_solution(solution, indices, wall_indices, bodies);
            
            u_avg(Nx*(i-1) + j,:) = compute_average_velocity(solution, solution_local, ...
                xmin, xmax, ymin, ymax, 400, 400) / (hx*hy);
            
            p_avg(Nx*(i-1) + j) = compute_average_pressure(solution, solution_local, ...
                xmin, xmax, ymin, ymax, 400, 400) / (hx*hy);
            
            p_grad_avg(Nx*(i-1) + j,:) = compute_average_pressure_gradient(solution, ...
                xmin, xmax, ymin, ymax, 400, 400) / (hx*hy);
            
            u_grad_avg(Nx*(i-1) + j,:, :) = compute_average_velocity_gradient(solution, ...
                xmin, xmax, ymin, ymax, 400, 400) / (hx*hy);
    end
end
end

function solution_local = create_local_solution(solution, indices, wall_indices, bodies)

if ~isempty(bodies)
    npan = (wall_indices(1,2)-wall_indices(1,1)+1)/16;
    panel_breaks = [];
    
    for k = bodies
        panel_breaks = [panel_breaks; solution.problem.domain.panel_breaks((k-1)*(npan+1)+1:k*(npan+1))];
    end
    solution_local = solution;
    solution_local.problem.domain.z = solution.problem.domain.z(indices);
    solution_local.problem.domain.zp = solution.problem.domain.zp(indices);
    solution_local.problem.domain.zpp = solution.problem.domain.zpp(indices);
    solution_local.problem.domain.wazp = solution.problem.domain.wazp(indices);
    solution_local.problem.domain.quad_weights = solution.problem.domain.quad_weights(indices);
    solution_local.problem.domain.wall_indices = wall_indices;
    solution_local.q = solution.q(indices,:);
    solution_local.local_indices = indices;
    solution_local.problem.domain.panel_breaks = panel_breaks;
else
    % No obstacles in cell
    solution_local.problem.domain.z = [];
    solution_local.problem.domain.zp = [];
    solution_local.problem.periodic = 1;
end
end