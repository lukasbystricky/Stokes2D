function sigma_avg = compute_pipe_average_stress(solution,xmin,xmax,ymin,ymax,npan)
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
sigma_avg = [0 0; 0 0];

% compute integral over boundary
[sxx,sxy,syx,syy] = evaluate_stress_on_surface(solution, solution, 'fluid');
t1 = n1.*sxx + n2.*syx;
t2 = n1.*sxy + n2.*syy;
sigma_avg(1,:) = sigma_avg(1,:) + [sum(x.*t1.*wazp) sum(x.*t2.*wazp)];
sigma_avg(2,:) = sigma_avg(2,:) + [sum(y.*t1.*wazp) sum(y.*t2.*wazp)];

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

[sxx,sxy,syx,syy] = evaluate_stress(solution, x, y);
t1 = n1.*sxx + n2.*syx;
t2 = n1.*sxy + n2.*syy;
sigma_avg(1,:) = sigma_avg(1,:) + [sum(x.*t1.*wazp) sum(x.*t2.*wazp)];
sigma_avg(2,:) = sigma_avg(2,:) + [sum(y.*t1.*wazp) sum(y.*t2.*wazp)];

% right side, n = (1,0)
right_side{1} = @(T) vertical_line(T,xmax,ymin,ymax,-1);
right_domain = discretize_domain(right_side,npan);
z = right_domain.z;
x = real(z);
y = imag(z);
wazp = right_domain.wazp;
n1 = ones(size(y));
n2 = zeros(size(y));

[sxx,sxy,syx,syy] = evaluate_stress(solution, x, y);
t1 = n1.*sxx + n2.*syx;
t2 = n1.*sxy + n2.*syy;
sigma_avg(1,:) = sigma_avg(1,:) + [sum(x.*t1.*wazp) sum(x.*t2.*wazp)];
sigma_avg(2,:) = sigma_avg(2,:) + [sum(y.*t1.*wazp) sum(y.*t2.*wazp)];

V = (xmax - xmin)*(ymax - ymin);
sigma_avg = sigma_avg/V;
end