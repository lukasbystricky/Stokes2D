function [] = plot_domain(problem)

figure()
hold on

z = problem.domain.z;
zp = problem.domain.zp;

Lx = problem.Lx;
Ly = problem.Ly;

plot(z, '-o');

% plot normal vector, recalling n = -1i*zp / |zp|
quiver(real(z), imag(z), real(-1i*zp)./abs(zp), imag(-1i*zp)./abs(zp));

xlim([min(real(z)), max(real(z))]);
ylim([min(imag(z)), max(imag(z))]);

rectangle('Position',[-Lx/2, -Lx/2, Lx, Ly], 'LineStyle', '--')

axis equal
