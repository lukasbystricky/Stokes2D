% Tests the slip boundary condition on the periodic Pouseille flow problem
% with curved walls.
%   1. Solve problem with non-zero slip length alpha.
%   2. Read of Dirichlet boundary condition.
%   3. Solve new Dirichlet problem.
%   4. Compare solutions, which now should be equal.

% Current state:
%   eta = inf
%       - alpha=0, amp=0 works
%       - alpha nnz, amp=0 works
%       - alpha=0, amp nnz works
%       - alpha=amp=1e-5 works
%       - alpha=amp=1e-1 works
%   eta = 1
%       - alpha=amp=1e-1 works (converges to 1e-10)

close all;
clearvars;
clc;
set(groot,'defaultAxesTickLabelInterpreter','latex');

% parameters
alpha = 1e-1;
amp = 0;

% create input structure
input_params = default_input_params('slip_pipe_test', 1);
%input_params.gmres_tol = 1e-9;

% modify structure as needed
input_params.box_size = [2,2];
input_params.h = 0.5;
input_params.panels = 10;
input_params.eta = inf;
input_params.plot_domain = 0;
input_params.slip = 1;

input_params.alpha = alpha;
input_params.A = amp;
input_params.d = 0.5;

problem_slip = sine_pipe_periodic(input_params);
input_params.alpha = 0;
problem_dirichlet = sine_pipe_periodic(input_params);

%% solve slip problem with curved boundary
N = length(problem_slip.domain.z);
rhs_slip = [zeros(2*N,1); problem_slip.pressure_gradient_x; problem_slip.pressure_gradient_y];
solution_slip = solve_stokes(problem_slip,rhs_slip,'fmm',0);
solution_slip.local_indices = 1:length(solution_slip.q);

%% read of boundary condition
domain_slip = problem_slip.domain;

% normalized tangent
nu1 = real(domain_slip.zp)./abs(domain_slip.zp);
nu2 = imag(domain_slip.zp)./abs(domain_slip.zp);

% normalized normal
n1 = real(-1i*domain_slip.zp)./abs(domain_slip.zp);
n2 = imag(-1i*domain_slip.zp)./abs(domain_slip.zp);

% on-surface velocity and stress
[u1_slip,u2_slip] = evaluate_velocity_on_surface(solution_slip,solution_slip);

% boundary condition in normal direction
udotn = u1_slip.*n1 + u2_slip.*n2;

% boundary condition in tangential direction
udotnu = u1_slip.*nu1 + u2_slip.*nu2;

%% solve new Dirichlet problem
g1 = udotn;
g2 = udotnu;
rhs_dirichlet = [g2; g1; problem_dirichlet.pressure_gradient_x; problem_dirichlet.pressure_gradient_y];
solution_dirichlet = solve_stokes(problem_dirichlet,rhs_dirichlet,'fmm',0);
solution_dirichlet.local_indices = 1:length(solution_dirichlet.q);

%% compute velocities
[Uslip,Vslip,X,Y] = evaluate_velocity(solution_slip,200,'fmm',0,'verbose',0);
[Udirichlet,Vdirichlet] = evaluate_velocity(solution_dirichlet,X,Y,'fmm',0,'verbose',0);
[u1_dirichlet,u2_dirichlet] = evaluate_velocity_on_surface(solution_dirichlet,solution_dirichlet);

%% plot off-surface velocities
figure('DefaultAxesFontSize',16);
subplot(2,2,1);
contourf(X,Y,Uslip,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$u$ slip','interpreter','latex','fontsize',16);

subplot(2,2,2);
contourf(X,Y,Vslip,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$v$ slip','interpreter','latex','fontsize',16);

subplot(2,2,3);
contourf(X,Y,Udirichlet,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$u$ Dirichlet','interpreter','latex','fontsize',16);

subplot(2,2,4);
contourf(X,Y,Vdirichlet,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$v$ Dirichlet','interpreter','latex','fontsize',16);
%saveas(gca,'/home/david/Documents/Research/kernels/slip_tests/pipe/vel_amp_1e-6_alpha_1e-5_eta_inf_sq_matlab','epsc');

%% plot on-surface velocities
figure('DefaultAxesFontSize',16);
subplot(1,2,1);
plot(g1);
xlabel('discretization points','interpreter','latex','fontsize',16);
ylabel('u1 slip','interpreter','latex','fontsize',16);
xlim([1,N]);
grid on;

subplot(1,2,2);
plot(g2);
xlabel('discretization points','interpreter','latex','fontsize',16);
ylabel('u2 slip','interpreter','latex','fontsize',16);
xlim([1,N]);
grid on;

%% plot off-surface errors
figure('DefaultAxesFontSize',16);
subplot(1,2,1);
contourf(X,Y,log10(abs(Uslip-Udirichlet)+eps),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
%caxis([-16,-12]);
axis equal;
title('$u$ Dirichlet vs slip: $\log_{10}$(absolute error)','interpreter','latex','fontsize',16);

subplot(1,2,2);
contourf(X,Y,log10(abs(Vslip-Vdirichlet)+eps),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
%caxis([-16,-12]);
axis equal;
title('$v$ Dirichlet vs slip: $\log_{10}$(absolute error)','interpreter','latex','fontsize',16);
%saveas(gca,'/home/david/Documents/Research/kernels/slip_tests/pipe/vel_err_amp_1e-1_alpha_1e-1_eta_inf_sq_matlab','epsc');

%% plot on-surface errors
figure('DefaultAxesFontSize',16);
subplot(2,1,1);
semilogy(abs(u1_slip-u1_dirichlet)+eps,'x--');
xlabel('discretization points','interpreter','latex','fontsize',16);
ylabel('absolute error','interpreter','latex','fontsize',16);
xlim([1,N]);
grid on;
title('$u$ on-surface Dirichlet vs slip','interpreter','latex','fontsize',16);

subplot(2,1,2);
semilogy(abs(u2_slip-u2_dirichlet)+eps,'x--');
xlabel('discretization points','interpreter','latex','fontsize',16);
ylabel('absolute error','interpreter','latex','fontsize',16);
xlim([1,N]);
grid on;
title('$v$ on-surface Dirichlet vs slip','interpreter','latex','fontsize',16);
%saveas(gca,'/home/david/Documents/Research/kernels/slip_tests/pipe/vel_err_amp_1e-1_alpha_1e-1_eta_inf_sq_matlab','epsc');
