function p_avg = compute_average_pressure(solution, solution_local, xmin, xmax, ...
    ymin, ymax, Nx, Ny)

% find quadrature points inside averaging volume
domain = solution_local.problem.domain;
z = domain.z;
zp = domain.zp;


p_avg = 0;
if ~isempty(z)
    
    x = real(z);
    y = imag(z);
    r = abs(z);
    n1 = real(-1i*zp)./abs(zp);
    n2 = imag(-1i*zp)./abs(zp);
    wazp = solution_local.problem.domain.wazp;
    
    %outward pointing normal
    n1 = -n1;
    n2 = -n2;
    
    % compute integral over boundary
    [Px, Py] = evaluate_pressure_gradient_on_surface(solution, solution_local);
    P = evaluate_pressure_on_surface(solution, solution_local);
    
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

V = (xmax - xmin)*(ymax - ymin);
p_avg = 0.25*p_avg/V;