function p_grad_avg = compute_pipe_average_pressure_gradient(solution,xmin,xmax,ymin,ymax,npan)
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
P = evaluate_pressure_on_surface(solution, solution,'fluid');

p_grad_avg = [0; 0];
p_grad_avg(1) = p_grad_avg(1) + sum(P.*n1.*wazp);
p_grad_avg(2) = p_grad_avg(2) + sum(P.*n2.*wazp);
%p_grad_avg = p_grad_avg + sum(P(1:end/2).*wazp(1:end/2))*[0;1] + sum(P(end/2+1:end).*wazp(end/2+1:end))*[0;-1];

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

P = evaluate_pressure(solution, x, y);
p_grad_avg(1) = p_grad_avg(1) + sum(P.*n1.*wazp);
p_grad_avg(2) = p_grad_avg(2) + sum(P.*n2.*wazp);
%p_grad_avg = p_grad_avg + sum(P.*wazp)*[-1;0];

% right side, n = (1,0)
right_side{1} = @(T) vertical_line(T,xmax,ymin,ymax,-1);
right_domain = discretize_domain(right_side,npan);
z = right_domain.z;
x = real(z);
y = imag(z);
wazp = right_domain.wazp;
n1 = ones(size(y));
n2 = zeros(size(y));

P = evaluate_pressure(solution, x, y);
p_grad_avg(1) = p_grad_avg(1) + sum(P.*n1.*wazp);
p_grad_avg(2) = p_grad_avg(2) + sum(P.*n2.*wazp);
%p_grad_avg = p_grad_avg + sum(P.*wazp)*[1;0];

V = (xmax - xmin)*(ymax - ymin);
p_grad_avg = p_grad_avg/V;
end