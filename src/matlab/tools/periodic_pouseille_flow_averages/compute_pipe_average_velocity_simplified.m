function u_avg = compute_pipe_average_velocity_simplified(solution,xmin,xmax,ymin,ymax,npan)
% quadrature points inside averaging volume
domain = solution.problem.domain;
z = domain.z;
zp = domain.zp;

% top and bottom boundary
x = real(z);
y = imag(z);
n1 = real(-1i*zp)./abs(zp);
n2 = imag(-1i*zp)./abs(zp);
wazp = domain.wazp;

% outward pointing normal
n1 = -n1;
n2 = -n2;

% initialize average velocity
u_avg = [0;0];

% compute integral over boundary
[U, V] = evaluate_velocity_on_surface(solution, solution);
u_avg(1) = u_avg(1) + sum((x.*U.*n1 + x.*V.*n2).*wazp);
u_avg(2) = u_avg(2) + sum((y.*V.*n1 + y.*U.*n2).*wazp);

% compute integral over the rest of the box using GL quadrature
solution.trim = 0;

% left side, n = (-1,0)
left_side{1} = @(T) vertical_line(T,xmin,ymin,ymax,1);
left_domain = discretize_domain(left_side,npan);
z = left_domain.z;
x = real(z);
y = imag(z);
wazp = left_domain.wazp;
n1 = -ones(size(y));
n2 = zeros(size(y));

[U, V] = evaluate_velocity(solution, x, y);
u_avg(1) = u_avg(1) + sum((x.*U.*n1 + x.*V.*n2).*wazp);
u_avg(2) = u_avg(2) + sum((y.*V.*n1 + y.*U.*n2).*wazp);

% right side, n = (1,0)
right_side{1} = @(T) vertical_line(T,xmax,ymin,ymax,-1);
right_domain = discretize_domain(right_side,npan);
z = right_domain.z;
x = real(z);
y = imag(z);
wazp = right_domain.wazp;
n1 = ones(size(y));
n2 = zeros(size(y));

[U, V] = evaluate_velocity(solution, x, y);
u_avg(1) = u_avg(1) + sum((x.*U.*n1 + x.*V.*n2).*wazp);
u_avg(2) = u_avg(2) + sum((y.*V.*n1 + y.*U.*n2).*wazp);

V = (xmax - xmin)*(ymax - ymin);
u_avg = u_avg/V;
end