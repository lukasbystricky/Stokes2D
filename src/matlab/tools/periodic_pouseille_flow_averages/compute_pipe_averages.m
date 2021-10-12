function [u_avg, u_grad_avg, p_avg, p_grad_avg] = compute_pipe_averages(solution, varargin)
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
u_avg = compute_pipe_average_velocity(solution,xmin,xmax,ymin,ymax,npan);
u_grad_avg = compute_pipe_average_velocity_gradient(solution,xmin,xmax,ymin,ymax,npan);
p_avg = compute_pipe_average_pressure(solution,xmin,xmax,ymin,ymax,npan);
p_grad_avg = compute_pipe_average_pressure_gradient(solution,xmin,xmax,ymin,ymax,npan);
end