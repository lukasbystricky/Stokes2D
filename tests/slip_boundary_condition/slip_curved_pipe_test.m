% Tests the slip boundary condition on the periodic Pouseille flow problem
% with curved boundaries.
%   1. Solve problem with non-zero slip length alpha.
%   2. Read of Dirichlet boundary condition.
%   3. Solve new Dirichlet problem.
%   4. Compare solutions, which now should be equal.

% Current state:
% GMRES does not converge for large alpha.

close all;
clearvars;
clc;
set(groot,'defaultAxesTickLabelInterpreter','latex');

% parameters
alpha = 1e-5;
amp = 1e-1;

% create input structure
input_params = default_input_params('slip_pipe_test', 1);

% modify structure as needed
input_params.box_size = [2,5];
input_params.h = 0.5;
input_params.panels = 20;
input_params.eta = inf;
input_params.plot_domain = 0;
input_params.slip = 1;

input_params.alpha = alpha;
input_params.A = amp;
input_params.d = 0.5;

problem_slip = sine_pipe_periodic(input_params);
input_params.alpha = 0;
problem_dirichlet = flat_pipe_periodic(input_params);

%% solve slip problem with curved boundary
N = length(problem_slip.domain.z);
rhs_slip = [zeros(2*N,1); problem_slip.pressure_gradient_x; problem_slip.pressure_gradient_y];
solution_slip = solve_stokes(problem_slip,rhs_slip,'fmm',0);

%% read of boundary condition
domain_slip = problem_slip.domain;

% normalized tangent
nu1 = real(domain_slip.zp)./abs(domain_slip.zp);
nu2 = imag(domain_slip.zp)./abs(domain_slip.zp);

% normalized normal
n1 = real(-1i*domain_slip.zp)./abs(domain_slip.zp);
n2 = imag(-1i*domain_slip.zp)./abs(domain_slip.zp);

% on-surface velocity and stress
solution_slip.local_indices = 1:N;
[u1_slip,u2_slip] = evaluate_velocity_on_surface(solution_slip, solution_slip);

% boundary condition in normal direction
udotn = u1_slip.*n1 + u2_slip.*n2;

% boundary condition in tangential direction
udotnu = u1_slip.*nu1 + u2_slip.*nu2;

%% solve new Dirichlet problem
g1 = udotn;
g2 = udotnu;
rhs_dirichlet = [g2; g1; problem_dirichlet.pressure_gradient_x; problem_dirichlet.pressure_gradient_y];
solution_dirichlet = solve_stokes(problem_dirichlet,rhs_dirichlet,'fmm',0);

%% compute velocities
[Uslip,Vslip,X,Y] = evaluate_velocity(solution_slip,200,'fmm',0,'verbose',0);
[Udirichlet,Vdirichlet] = evaluate_velocity(solution_slip,X,Y,'fmm',0,'verbose',0);

%% plot solutions
figure('DefaultAxesFontSize',16);
subplot(2,2,1);
contourf(X,Y,Udirichlet,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
caxis('auto');
lim = caxis;
axis equal;
title('$u$ Dirichlet','interpreter','latex','fontsize',16);

subplot(2,2,2);
contourf(X,Y,Vdirichlet,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$v$ Dirichlet','interpreter','latex','fontsize',16);

subplot(2,2,3);
contourf(X,Y,Uslip,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
caxis(lim);
axis equal;
title('$u$ slip','interpreter','latex','fontsize',16);

subplot(2,2,4);
contourf(X,Y,Vslip,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$v$ slip','interpreter','latex','fontsize',16);
%saveas(gca,'/home/david/Documents/Research/kernels/slip_tests/pipe/vel_amp_1e-6_alpha_1e-5_eta_inf_sq_matlab','epsc');

%% plot errors
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
