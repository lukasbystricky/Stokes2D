% Taylor-Couette problem with slip boundary condition. Analytical solution
% is with slip length -alpha by our definition. Therefore, we have an extra
% minus sign in the analytical solution.

close all;
clearvars;
clc;
set(groot,'defaultAxesTickLabelInterpreter','latex');
format long;

% create input structure
input_params = default_input_params('couette', 0);

% modify structure as needed, or add additional problem-dependent params
input_params.panels = [10, 10];
input_params.plot_domain = 0;

% set up radii and centers, centers given as complex numbers (x+iy)
input_params.radii = [1, 0.5]; % the outer wall must be the first radii
input_params.centers = [0, 0];
input_params.omega = 1; % angular velocity of inner wall
input_params.slip = 1;
input_params.alpha = 1e-5;

% setup problem
problem = couette(input_params);

% quantities for exact solution
ro = input_params.radii(1);
ri = input_params.radii(2);
omega = input_params.omega;
alpha = -input_params.alpha;    % note minus sign

% exact solutions
A = (1-2*alpha/ro)/ro^2;
B = (1+2*alpha/ri)/ri^2;
exact_solution_r = @(x,y) zeros(size(x));
exact_solution_theta = @(x,y) (omega/(A-B)) * (A*sqrt(x.^2 + y.^2) - 1./sqrt(x.^2 + y.^2));

%% solve the problem
z = problem.domain.z;
N = length(z);

% setup rhs
nu1 = real(problem.domain.zp)./abs(problem.domain.zp);
nu2 = imag(problem.domain.zp)./abs(problem.domain.zp);
n1 = real(-1i*problem.domain.zp)./abs(problem.domain.zp);
n2 = imag(-1i*problem.domain.zp)./abs(problem.domain.zp);
% rhs = [-alpha*(2*omega/(A-B)/ro^2)*ones(N/2,1);
%        ri*omega*ones(N/2,1)+alpha*(2*omega/(A-B)/ri^2)*ones(N/2,1);
%        zeros(N,1)];
rhs = [zeros(N/2,1);
       ri*omega*ones(N/2,1);
       zeros(N,1)];

% solve the problem
solution = solve_stokes(problem,rhs,'fmm',0);
solution.local_indices = 1:length(solution.q);

%% compute velocity and stress
% grid in polar coordinates with M number of radial and angular points
M = 200;
r = linspace(input_params.radii(2)+1e-2, input_params.radii(1) - 1e-2, M);
h = (2*pi)/M;
theta = (0:h:2*pi);
[R, T] = meshgrid(r, theta);
X = R.*cos(T);
Y = R.*sin(T);
x = real(problem.domain.z);
y = imag(problem.domain.z);

[U,V,X,Y] = evaluate_velocity(solution,X,Y,'fmm',0,'verbose',0);
[usurf,vsurf] = evaluate_velocity_on_surface(solution,solution);
[sxxsurf,sxysurf,syxsurf,syysurf] = evaluate_stress_on_surface(solution,solution,'fluid');

% convert Cartesian velocity to radial and angular velocity, using 
% relationships:
% theta_hat = -i sin(theta) +j cos(theta) = i y/r + j x/r
% r_hat = i cos(theta) + j sin(theta) = i x/r + j y/r
Ur = (U.*X + V.*Y)./sqrt(X.^2 + Y.^2);
Utheta = (-U.*Y + V.*X)./sqrt(X.^2 + Y.^2);
ursurf = (usurf.*x + vsurf.*y)./sqrt(x.^2 + y.^2);
uthetasurf = (-usurf.*y + vsurf.*x)./sqrt(x.^2 + y.^2);

% check boundary values
udotn = usurf.*n1 + vsurf.*n2;
udotnu = usurf.*nu1 + vsurf.*nu2;
t1 = n1.*sxxsurf + n2.*syxsurf;
t2 = n1.*sxysurf + n2.*syysurf;
nudott = t1.*nu1 + t2.*nu2;

b1 = udotn;
b2 = udotnu - alpha*nudott;

%% plot solutions
figure('DefaultAxesFontSize',16);
subplot(2,2,1);
contourf(X,Y,Utheta,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$u_\theta$','interpreter','latex','fontsize',16);

subplot(2,2,2);
contourf(X,Y,exact_solution_theta(X,Y),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('exact $u_\theta$','interpreter','latex','fontsize',16);

subplot(2,2,3);
contourf(X,Y,Ur,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$u_r$','interpreter','latex','fontsize',16);

%% plot errors
figure('DefaultAxesFontSize',16);
subplot(1,2,1);
contourf(X,Y,log10(abs(Ur)+eps),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
%caxis([-16,-12]);
axis equal;
title('$u_r$: $\log_{10}$(absolute error)','interpreter','latex','fontsize',16);

subplot(1,2,2);
contourf(X,Y,log10(abs((Utheta-exact_solution_theta(X,Y))./max(max(abs(exact_solution_theta(X,Y)))))+eps),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
%caxis([-16,-12]);
axis equal;
title('$u_\theta$: $\log_{10}$(relative error)','interpreter','latex','fontsize',16);
%saveas(gca,'/home/david/Documents/Research/kernels/slip_tests/couette/domain_u_err','epsc');

figure('DefaultAxesFontSize',16);
semilogy(abs(ursurf));
title('$u_r$ on-surface: abs error','interpreter','latex','fontsize',16);
grid on;
xlim([1,length(x)]);

figure('DefaultAxesFontSize',16);
semilogy(abs(b1));
title('udotn on-surface: abs error','interpreter','latex','fontsize',16);
grid on;
xlim([1,length(x)]);
%saveas(gca,'/home/david/Documents/Research/kernels/slip_tests/couette/u_surf_err','epsc');

figure('DefaultAxesFontSize',16);
semilogy([abs(b2(1:end/2)-0); abs(b2(end/2+1:end)-0.5)]);
title('udotnu + alpha*nudott on-surface: abs error','interpreter','latex','fontsize',16);
grid on;
xlim([1,length(x)]);
