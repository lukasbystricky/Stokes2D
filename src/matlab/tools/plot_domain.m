function [] = plot_domain(problem)

figure()
hold on

Lx = problem.Lx;
Ly = problem.Ly;

x = mod(real(problem.domain.z)+Lx/2,Lx)-Lx/2;
y = mod(imag(problem.domain.z)+Ly/2,Ly)-Ly/2;

z = x + 1i*y;
zp = problem.domain.zp;

plot(z, '-o');

% plot normal vector, recalling n = -1i*zp / |zp|
quiver(real(z), imag(z), real(-1i*zp)./abs(zp), imag(-1i*zp)./abs(zp));

xlim([min(real(z)), max(real(z))]);
ylim([min(imag(z)), max(imag(z))]);

rectangle('Position',[-Lx/2, -Ly/2, Lx, Ly], 'LineStyle', '--')

axis equal
