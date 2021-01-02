function u_grad_avg = compute_average_velocity_gradient(solution, xmin, xmax, ...
    ymin, ymax, Nx, Ny, varargin)

u_grad_avg = [0 0; 0,0];

% compute integral over box boundary
hx = (xmax - xmin)/Nx;
hy = (ymax - ymin)/Ny;

% left side, n = (-1,0)
y = ymin : hy : ymax;
x = xmin*ones(size(y));
w = hy*ones(size(x));
w(1) = hy/2;
w(end) = hy/2;

[u1, u2] = evaluate_velocity(solution, x, y);
u_grad_avg = u_grad_avg + [-sum(w.*u1), -sum(w.*u2); 0, 0];

% right side, n = (1,0)
y = ymin : hy : ymax;
x = xmax*ones(size(y));
w = hy*ones(size(x));
w(1) = hy/2;
w(end) = hy/2;

[u1, u2] = evaluate_velocity(solution, x, y);
u_grad_avg = u_grad_avg + [sum(w.*u1), sum(w.*u2); 0, 0];

% top side, n = (0,1)
x = xmin : hx : xmax;
y = ymax*ones(size(x));
w = hx*ones(size(x));
w(1) = hx/2;
w(end) = hx/2;

[u1, u2] = evaluate_velocity(solution, x, y);
u_grad_avg = u_grad_avg + [0, 0; sum(w.*u1), sum(w.*u2)];

% bottom side, n = (0, -1)
x = xmin : hx : xmax;
y = ymin*ones(size(x));
w = hx*ones(size(x));
w(1) = hx/2;
w(end) = hx/2;

[u1, u2] = evaluate_velocity(solution, x, y);
u_grad_avg = u_grad_avg + [0, 0; -sum(w.*u1), -sum(w.*u2)];