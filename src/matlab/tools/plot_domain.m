function [] = plot_domain(problem, varargin)

if nargin == 1
    figure()
else
    figure(varargin{1});
end

hold on
z = problem.domain.z;
zp = problem.domain.zp;

if problem.periodic
    Lx = problem.Lx;
    Ly = problem.Ly;
    
    x = mod(real(problem.domain.z)+Lx/2,Lx)-Lx/2;
    y = mod(imag(problem.domain.z)+Ly/2,Ly)-Ly/2;
    z = x + 1i*y;
end

plot(z, '-o');

% plot normal vector, recalling n = -1i*zp / |zp|
quiver(real(z), imag(z), real(-1i*zp)./abs(zp), imag(-1i*zp)./abs(zp));

xlim([min(real(z)), max(real(z))]);
ylim([min(imag(z)), max(imag(z))]);

if problem.periodic
    rectangle('Position',[-Lx/2, -Ly/2, Lx, Ly], 'LineStyle', '--')
    xlim([-1.1*Lx/2, 1.1*Lx/2]);
    ylim([-1.1*Ly/2, 1.1*Ly/2]);
end

axis equal
