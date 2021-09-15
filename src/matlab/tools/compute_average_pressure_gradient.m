function p_grad_avg = compute_average_pressure_gradient(solution, xmin, xmax, ...
    ymin, ymax, Nx, Ny)

p_grad_avg = [0;0];

% compute integral over box boundary
hx = (xmax - xmin)/Nx;
hy = (ymax - ymin)/Ny;

% left side, n = (-1,0)
y = ymin : hy : ymax;
x = xmin*ones(size(y));
w = hy*ones(size(x));
w(1) = hy/2;
w(end) = hy/2;

P = evaluate_pressure(solution, x, y);
p_grad_avg = p_grad_avg + sum(P.*w)*[-1;0];

% right side, n = (1,0)
y = ymin : hy : ymax;
x = xmax*ones(size(y));
w = hy*ones(size(x));
w(1) = hy/2;
w(end) = hy/2;

P = evaluate_pressure(solution, x, y);
p_grad_avg = p_grad_avg + sum(P.*w)*[1;0];

% top side, n = (0,1)
x = xmin : hx : xmax;
y = ymax*ones(size(x));
w = hx*ones(size(x));
w(1) = hx/2;
w(end) = hx/2;

P = evaluate_pressure(solution, x, y);
p_grad_avg = p_grad_avg + sum(P.*w)*[0;1];

% bottom side, n = (0, -1)
x = xmin : hx : xmax;
y = ymin*ones(size(x));
w = hx*ones(size(x));
w(1) = hx/2;
w(end) = hx/2;

P = evaluate_pressure(solution, x, y);
p_grad_avg = p_grad_avg + sum(P.*w)*[0;-1];

V = (xmax - xmin)*(ymax - ymin);
p_grad_avg = p_grad_avg / V;