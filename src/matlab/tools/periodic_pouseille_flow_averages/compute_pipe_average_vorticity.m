function omega_avg = compute_pipe_average_vorticity(solution,xmin,xmax,ymin,ymax,npan)
% quadrature points inside averaging volume
domain = solution.problem.domain;
z = domain.z;
zp = domain.zp;

% top and bottom boundary
n1 = real(-1i*zp)./abs(zp);
n2 = imag(-1i*zp)./abs(zp);
wazp = domain.wazp;

% outward pointing normal
n1 = -n1;
n2 = -n2;

omega_avg = 0;

% compute integral over boundary
[U, V] = evaluate_velocity_on_surface(solution, solution);
omega_avg = omega_avg + sum((n1.*V-U.*n2).*wazp);

% compute integral over the rest of the box using GL quadrature
solution.trim = 0;

% left side, n = (-1,0)
left_side{1} = @(T) vertical_line(T,xmin,ymin,ymax,1);
left_domain = discretize_domain(left_side,npan);
z = left_domain.z;
x = real(z);
y = imag(z);
wazp = left_domain.wazp;
n1 = -ones(size(z));
n2 = zeros(size(z));

[U, V] = evaluate_velocity(solution, x, y);
omega_avg = omega_avg + sum((n1.*V-U.*n2).*wazp);

% right side, n = (1,0)
right_side{1} = @(T) vertical_line(T,xmax,ymin,ymax,-1);
right_domain = discretize_domain(right_side,npan);
z = left_domain.z;
x = real(z);
y = imag(z);
wazp = right_domain.wazp;
n1 = ones(size(z));
n2 = zeros(size(z));

[U, V] = evaluate_velocity(solution, x, y);
omega_avg = omega_avg + sum((n1.*V-U.*n2).*wazp);

V = (xmax - xmin)*(ymax - ymin);
omega_avg = omega_avg/V;
end