function u_avg = compute_average_velocity(solution, xmin, xmax, ...
    ymin, ymax, Nx, Ny, bodies)

% find quadrature points inside averaging volume
indices_bodies = [];
for i = bodies
    indices = [indices_bodies, solution.problem.domain.wall_indices(bodies,1) : ...
                solution.problem.domain.wall_indices(bodies,2)];
end

z = solution.problem.domain.z(indices);
zp = solution.problem.domain.zp(indices);

u_avg = [0;0];

if ~isempty(indices)
    
    x = real(z(indices));
    y = imag(z(indices));
    r = abs(z(indices));
    n1 = real(-1i*zp(indices))./abs(zp(indices));
    n2 = imag(-1i*zp(indices))./abs(zp(indices));
    wazp = solution.problem.domain.wazp(indices);
    
    %outward pointing normal
    n1 = -n1;
    n2 = -n2;
    
    
    % compute integral over boundary
    [Ux, Uy, Vx, Vy] = evaluate_velocity_gradient_on_surface(solution, bodies);
    [U, V] = evaluate_velocity_on_surface(solution, bodies);
    
    Ux = Ux(indices);
    Uy = Uy(indices);
    Vx = Vx(indices);
    Vy = Vy(indices);
    U = U(indices);
    V = V(indices);
    
    
    u_avg(1) = u_avg(1) + sum((r.^2.*(Ux.*n1 + Vx.*n2) - 2*(x.*U + y.*V).*n1).*wazp);
    u_avg(2) = u_avg(2) + sum((r.^2.*(Uy.*n1 + Vy.*n2) - 2*(x.*U + y.*V).*n2).*wazp);
end

% compute integral over box boundary
hx = (xmax - xmin)/Nx;
hy = (ymax - ymin)/Ny;

% left side, n = (-1,0)
y = ymin : hy : ymax;
x = xmin*ones(size(y));
r = sqrt(x.^2 + y.^2);
n1 = -ones(size(y));
n2 = zeros(size(y));
w = hy*ones(size(x));
w(1) = hy/2;
w(end) = hy/2;

[Ux, Uy, Vx, Vy] = evaluate_velocity_gradient(solution, x, y);
[U, V] = evaluate_velocity(solution, x, y);
u_avg(1) = u_avg(1) + sum((r.^2.*(Ux.*n1 + Vx.*n2) - 2*(x.*U + y.*V).*n1).*w);
u_avg(2) = u_avg(2) + sum((r.^2.*(Uy.*n1 + Vy.*n2) - 2*(x.*U + y.*V).*n2).*w);

% right side, n = (1,0)
y = ymin : hy : ymax;
x = xmax*ones(size(y));
r = sqrt(x.^2 + y.^2);
n1 = ones(size(y));
n2 = zeros(size(y));
w = hy*ones(size(x));
w(1) = hy/2;
w(end) = hy/2;

[Ux, Uy, Vx, Vy] = evaluate_velocity_gradient(solution, x, y);
[U, V] = evaluate_velocity(solution, x, y);
u_avg(1) = u_avg(1) + sum((r.^2.*(Ux.*n1 + Vx.*n2) - 2*(x.*U + y.*V).*n1).*w);
u_avg(2) = u_avg(2) + sum((r.^2.*(Uy.*n1 + Vy.*n2) - 2*(x.*U + y.*V).*n2).*w);

% top side, n = (0,1)
x = xmin : hx : xmax;
y = ymax*ones(size(x));
r = sqrt(x.^2 + y.^2);
n1 = zeros(size(y));
n2 = ones(size(y));
w = hx*ones(size(x));
w(1) = hx/2;
w(end) = hx/2;

[Ux, Uy, Vx, Vy] = evaluate_velocity_gradient(solution, x, y);
[U, V] = evaluate_velocity(solution, x, y);
u_avg(1) = u_avg(1) + sum((r.^2.*(Ux.*n1 + Vx.*n2) - 2*(x.*U + y.*V).*n1).*w);
u_avg(2) = u_avg(2) + sum((r.^2.*(Uy.*n1 + Vy.*n2) - 2*(x.*U + y.*V).*n2).*w);

% bottom side, n = (0, -1)
x = xmin : hx : xmax;
y = ymin*ones(size(x));
r = sqrt(x.^2 + y.^2);
n1 = zeros(size(y));
n2 = -ones(size(y));
w = hx*ones(size(x));
w(1) = hx/2;
w(end) = hx/2;

[Ux, Uy, Vx, Vy] = evaluate_velocity_gradient(solution, x, y);
[U, V] = evaluate_velocity(solution, x, y);
u_avg(1) = u_avg(1) + sum((r.^2.*(Ux.*n1 + Vx.*n2) - 2*(x.*U + y.*V).*n1).*w);
u_avg(2) = u_avg(2) + sum((r.^2.*(Uy.*n1 + Vy.*n2) - 2*(x.*U + y.*V).*n2).*w);

u_avg = -0.5*u_avg;