function p_avg = compute_pipe_average_pressure(solution,xmin,xmax,ymin,ymax,npan)
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
[Px, Py] = evaluate_pressure_gradient_on_surface(solution, solution,'fluid');
P = evaluate_pressure_on_surface(solution, solution,'fluid');

p_avg = 0;
p_avg = p_avg + sum((2*P.*(x.*n1+y.*n2) - r.^2.*(Px.*n1 + Py.*n2)).*wazp);

% compute integral over the rest of the box using GL quadrature
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

[Px, Py] = evaluate_pressure_gradient(solution,x,y);
P = evaluate_pressure(solution,x,y);
p_avg = p_avg + sum((2*P.*(x.*n1+y.*n2) - r.^2.*(Px.*n1 + Py.*n2)).*wazp);

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

[Px, Py] = evaluate_pressure_gradient(solution,x,y);
P = evaluate_pressure(solution,x,y);
p_avg = p_avg + sum((2*P.*(x.*n1+y.*n2) - r.^2.*(Px.*n1 + Py.*n2)).*wazp);

V = (xmax - xmin)*(ymax - ymin);
p_avg = 0.25*p_avg/V;
end