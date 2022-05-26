function [u_avg_alt1, u_avg_alt2, u_avg_alt3, u_grad_avg, p_avg, p_grad_avg, omega_avg, stress_avg] = compute_pipe_averages(solution, varargin)
% number of GL panels on left and right side of box
if nargin > 1
    npan = varargin{1};
else
    npan = 1;
end

% averaging volume limits
xmin = -solution.problem.Lx/2;
xmax = solution.problem.Lx/2;
ymin = -solution.problem.h;
ymax = solution.problem.h;

% needed for on-surface evaluation
solution.local_indices = 1:length(solution.q);

% compute averages (\int_V X dV)/V
u_avg_alt1 = compute_pipe_average_velocity(solution,xmin,xmax,ymin,ymax,npan);
u_avg_alt2 = compute_pipe_average_velocity_using_stress(solution,xmin,xmax,ymin,ymax,npan);
u_avg_alt3 = compute_pipe_average_velocity_simplified(solution,xmin,xmax,ymin,ymax,npan);
u_grad_avg = compute_pipe_average_velocity_gradient(solution,xmin,xmax,ymin,ymax,npan);
p_avg = compute_pipe_average_pressure(solution,xmin,xmax,ymin,ymax,npan);
p_grad_avg = compute_pipe_average_pressure_gradient(solution,xmin,xmax,ymin,ymax,npan);
omega_avg = compute_pipe_average_vorticity(solution,xmin,xmax,ymin,ymax,npan);
stress_avg = compute_pipe_average_stress(solution,xmin,xmax,ymin,ymax,npan);
end