function p_grad_avg = compute_average_pressure_gradient(solution, xmin, xmax, ...
    ymin, ymax, Nx, Ny, varargin)

p_grad_avg = [0;0];

% find quadrature points inside averaging volume
z = solution.problem.domain.z;
zp = solution.problem.domain.zp;
indices = 1:length(z);
indices = indices(real(z) > xmin & real(z) < xmax & ...
    imag(z) > ymin & imag(z) < ymax);

% if ~isempty(indices)
%     
%     x = real(z(indices));
%     y = imag(z(indices));
%     r = abs(z(indices));
%     n1 = real(-1i*zp(indices))./abs(zp(indices));
%     n2 = imag(-1i*zp(indices))./abs(zp(indices));
%     wazp = solution.problem.domain.wazp(indices);
%     
%     %outer normal
%     n1 = -n1;
%     n2 = -n2;
%     
%     if nargin == 7   
%         
%         % compute integral over boundary
%         P = evaluate_pressure_on_surface(solution);
%         
%         P = P(indices);       
%     else
%         
%         % integrals over boundary can also be passed in
%         P = varargin{1}(indices);     
%         
%     end
%     
%     p_grad_avg(1) = p_grad_avg(1) + sum(P.*n1.*wazp);
%     p_grad_avg(2) = p_grad_avg(2) + sum(P.*n2.*wazp);
% end

% compute integral over box boundary
hx = (xmax - xmin)/Nx;
hy = (ymax - ymin)/Ny;

% left side, n = (-1,0)
y = ymin : hy : ymax;
x = xmin*ones(size(y));
w = hy*ones(size(x));
w(1) = hy/2;
w(end) = hy/2;

P = evaluate_pressure(solution, x, y);
p_grad_avg = p_grad_avg + sum(P.*w)*[-1;0];

% right side, n = (1,0)
y = ymin : hy : ymax;
x = xmax*ones(size(y));
w = hy*ones(size(x));
w(1) = hy/2;
w(end) = hy/2;

P = evaluate_pressure(solution, x, y);
p_grad_avg = p_grad_avg + sum(P.*w)*[1;0];

% top side, n = (0,1)
x = xmin : hx : xmax;
y = ymax*ones(size(x));
w = hx*ones(size(x));
w(1) = hx/2;
w(end) = hx/2;

P = evaluate_pressure(solution, x, y);
p_grad_avg = p_grad_avg + sum(P.*w)*[0;1];

% bottom side, n = (0, -1)
x = xmin : hx : xmax;
y = ymin*ones(size(x));
w = hx*ones(size(x));
w(1) = hx/2;
w(end) = hx/2;

P = evaluate_pressure(solution, x, y);
p_grad_avg = p_grad_avg + sum(P.*w)*[0;-1];