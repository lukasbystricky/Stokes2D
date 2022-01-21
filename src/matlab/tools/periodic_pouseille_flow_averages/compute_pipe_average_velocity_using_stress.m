function u_avg = compute_pipe_average_velocity_using_stress(solution,xmin,xmax,ymin,ymax,npan)
% currently not working, probably something wrong in the volume integral
% formula

% quadrature points inside averaging volume
domain = solution.problem.domain;
z = domain.z;
zp = domain.zp;

% top and bottom boundary
x = real(z);
y = imag(z);
r = abs(z);
n1 = real(-1i*zp)./abs(zp);
n2 = imag(-1i*zp)./abs(zp);
wazp = domain.wazp;

% outward pointing normal
n1 = -n1;
n2 = -n2;

% compute integral over boundary
[U, V] = evaluate_velocity_on_surface(solution, solution);
P = evaluate_pressure_on_surface(solution, solution, 'fluid');
[Px, Py] = evaluate_pressure_gradient_on_surface(solution, solution, 'fluid');
[Sxx, Sxy, Syx, Syy] = evaluate_stress_on_surface(solution, solution, 'fluid');
T1 = n1.*Sxx + n2.*Sxy;
T2 = n1.*Syx + n2.*Syy;

u_avg = [0;0];
u_avg(1) = u_avg(1) + sum(((2*(U.*(x.*n1 + y.*n2) + n1.*(U.*x + V.*y)) - r.^2.*T1)/6 - ...
    (2*P.*x.*(x.*n1 + y.*n2) + 2*r.^2.*P.*n1 - r.^2.*x.*(Px.*n1 + Py.*n2))/24).*wazp);
u_avg(2) = u_avg(2) + sum(((2*(V.*(x.*n1 + y.*n2) + n2.*(U.*x + V.*y)) - r.^2.*T2)/6 - ...
    (2*P.*y.*(x.*n1 + y.*n2) + 2*r.^2.*P.*n2 - r.^2.*y.*(Px.*n1 + Py.*n2))/24).*wazp);

% compute integral over the rest of the box using GL quadrature
solution.trim = 0;

% left side, n = (-1,0)
left_side{1} = @(T) vertical_line(T,xmin,ymin,ymax,1);
left_domain = discretize_domain(left_side,npan);
z = left_domain.z;
x = real(z);
y = imag(z);
r = abs(z);
wazp = left_domain.wazp;
n1 = -ones(size(y));
n2 = zeros(size(y));

[U, V] = evaluate_velocity(solution, x, y);
P = evaluate_pressure(solution, x, y);
[Px, Py] = evaluate_pressure_gradient(solution, x, y);
[Sxx, Sxy, Syx, Syy] = evaluate_stress(solution, x, y);
T1 = n1.*Sxx + n2.*Sxy;
T2 = n1.*Syx + n2.*Syy;
u_avg(1) = u_avg(1) + sum(((2*(U.*(x.*n1 + y.*n2) + n1.*(U.*x + V.*y)) - r.^2.*T1)/6 - ...
    (2*P.*x.*(x.*n1 + y.*n2) + 2*r.^2.*P.*n1 - r.^2.*x.*(Px.*n1 + Py.*n2))/24).*wazp);
u_avg(2) = u_avg(2) + sum(((2*(V.*(x.*n1 + y.*n2) + n2.*(U.*x + V.*y)) - r.^2.*T2)/6 - ...
    (2*P.*y.*(x.*n1 + y.*n2) + 2*r.^2.*P.*n2 - r.^2.*y.*(Px.*n1 + Py.*n2))/24).*wazp);

% right side, n = (1,0)
right_side{1} = @(T) vertical_line(T,xmax,ymin,ymax,-1);
right_domain = discretize_domain(right_side,npan);
z = right_domain.z;
x = real(z);
y = imag(z);
r = abs(z);
wazp = right_domain.wazp;
n1 = ones(size(y));
n2 = zeros(size(y));

[U, V] = evaluate_velocity(solution, x, y);
P = evaluate_pressure(solution, x, y);
[Px, Py] = evaluate_pressure_gradient(solution, x, y);
[Sxx, Sxy, Syx, Syy] = evaluate_stress(solution, x, y);
T1 = n1.*Sxx + n2.*Sxy;
T2 = n1.*Syx + n2.*Syy;
u_avg(1) = u_avg(1) + sum(((2*(U.*(x.*n1 + y.*n2) + n1.*(U.*x + V.*y)) - r.^2.*T1)/6 - ...
    (2*P.*x.*(x.*n1 + y.*n2) + 2*r.^2.*P.*n1 - r.^2.*x.*(Px.*n1 + Py.*n2))/24).*wazp);
u_avg(2) = u_avg(2) + sum(((2*(V.*(x.*n1 + y.*n2) + n2.*(U.*x + V.*y)) - r.^2.*T2)/6 - ...
    (2*P.*y.*(x.*n1 + y.*n2) + 2*r.^2.*P.*n2 - r.^2.*y.*(Px.*n1 + Py.*n2))/24).*wazp);

V = (xmax - xmin)*(ymax - ymin);
u_avg = u_avg/V;
end