function u_grad_avg = compute_pipe_average_velocity_gradient(solution,xmin,xmax,ymin,ymax,npan)
% quadrature points inside averaging volume
domain = solution.problem.domain;
zp = domain.zp;

% top and bottom boundary
n1 = real(-1i*zp)./abs(zp);
n2 = imag(-1i*zp)./abs(zp);
wazp = domain.wazp;

% outward pointing normal
n1 = -n1;
n2 = -n2;

% compute integral over boundary
[U, V] = evaluate_velocity_on_surface(solution, solution);

u_grad_avg = [0 0; 0 0];
u_grad_avg(1,:) = u_grad_avg(1,:) + [sum(n1.*U.*wazp), sum(n1.*V.*wazp)];
u_grad_avg(2,:) = u_grad_avg(2,:) + [sum(n2.*U.*wazp), sum(n2.*V.*wazp)];

% compute integral over the rest of the box using GL quadrature
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
u_grad_avg(1,:) = u_grad_avg(1,:) + [sum(n1.*U.*wazp), sum(n1.*V.*wazp)];
u_grad_avg(2,:) = u_grad_avg(2,:) + [sum(n2.*U.*wazp), sum(n2.*V.*wazp)];

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
u_grad_avg(1,:) = u_grad_avg(1,:) + [sum(n1.*U.*wazp), sum(n1.*V.*wazp)];
u_grad_avg(2,:) = u_grad_avg(2,:) + [sum(n2.*U.*wazp), sum(n2.*V.*wazp)];

V = (xmax - xmin)*(ymax - ymin);
u_grad_avg = u_grad_avg/V;
end