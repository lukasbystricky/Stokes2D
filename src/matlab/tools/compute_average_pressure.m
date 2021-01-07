function p_avg = compute_average_pressure(solution, xmin, xmax, ...
    ymin, ymax, Nx, Ny, bodies)

% find quadrature points inside averaging volume

indices_bodies = [];
for i = bodies
    indices = [indices_bodies, solution.problem.domain.wall_indices(bodies,1) : ...
                solution.problem.domain.wall_indices(bodies,2)];
end

z = solution.problem.domain.z(indices);
zp = solution.problem.domain.zp(indices);


p_avg = 0;
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
    [Px, Py] = evaluate_pressure_gradient_on_surface(solution, bodies);
    P = evaluate_pressure_on_surface(solution, bodies);
    
    Px = Px(indices);
    Py = Py(indices);
    P = P(indices);
    
    
    
    p_avg = p_avg + sum((2*P.*(x.*n1+y.*n2) - r.^2.*(Px.*n1 + Py.*n2)).*wazp);
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

[Px, Py] = evaluate_pressure_gradient(solution,x,y);
P = evaluate_pressure(solution,x,y);
p_avg = p_avg + sum((2*P.*(x.*n1+y.*n2) - r.^2.*(Px.*n1 + Py.*n2)).*w);

% right side, n = (1,0)
y = ymin : hy : ymax;
x = xmax*ones(size(y));
r = sqrt(x.^2 + y.^2);
n1 = ones(size(y));
n2 = zeros(size(y));
w = hy*ones(size(x));
w(1) = hy/2;
w(end) = hy/2;

[Px, Py] = evaluate_pressure_gradient(solution,x,y);
P = evaluate_pressure(solution,x,y);
p_avg = p_avg + sum((2*P.*(x.*n1+y.*n2) - r.^2.*(Px.*n1 + Py.*n2)).*w);

% top side, n = (0,1)
x = xmin : hx : xmax;
y = ymax*ones(size(x));
r = sqrt(x.^2 + y.^2);
n1 = zeros(size(y));
n2 = ones(size(y));
w = hx*ones(size(x));
w(1) = hx/2;
w(end) = hx/2;

[Px, Py] = evaluate_pressure_gradient(solution,x,y);
P = evaluate_pressure(solution,x,y);
p_avg = p_avg + sum((2*P.*(x.*n1+y.*n2) - r.^2.*(Px.*n1 + Py.*n2)).*w);

% bottom side, n = (0, -1)
x = xmin : hx : xmax;
y = ymin*ones(size(x));
r = sqrt(x.^2 + y.^2);
n1 = zeros(size(y));
n2 = -ones(size(y));
w = hx*ones(size(x));
w(1) = hx/2;
w(end) = hx/2;

[Px, Py] = evaluate_pressure_gradient(solution,x,y);
P = evaluate_pressure(solution,x,y);
p_avg = p_avg + sum((2*P.*(x.*n1+y.*n2) - r.^2.*(Px.*n1 + Py.*n2)).*w);

p_avg = 0.25*p_avg;