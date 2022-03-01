% Test to compare the new matvec that can handle slip boundary condition
% with the old matvec that can only do Dirichlet. We do this by setting the
% slip length alpha as well as the corresponding boundary conditions to be
% zero. The flow is pressure driven and the geometry is a periodic curved
% pipe.

close all;
clearvars;
clc;
set(groot,'defaultAxesTickLabelInterpreter','latex');

% create input structure
input_params = default_input_params('no_slip_curved_pipe_matvec', 1);

% modify structure as needed
input_params.box_size = [2,5];
input_params.h = 0.5;
input_params.panels = 10;
input_params.eta = 1;
input_params.plot_domain = 0;
input_params.alpha = 0;
input_params.A = 1e-1;
input_params.d = 0.5;

problem = sine_pipe_periodic(input_params);

%% solve problem using new and old matvec with no slip
N = length(problem.domain.z);
rhs = [zeros(2*N,1); problem.pressure_gradient_x; problem.pressure_gradient_y];

problem.slip = 0;
solution_old = solve_stokes(problem,rhs,'fmm',0);
solution_old.local_indices = 1:length(solution_old.q);

problem.slip = 1;
solution_new = solve_stokes(problem,rhs,'fmm',0);
solution_new.local_indices = 1:length(solution_new.q);

%% compute velocities
[Uold,Vold,X,Y] = evaluate_velocity(solution_old,500,'fmm',1,'verbose',0);
[Unew,Vnew] = evaluate_velocity(solution_new,X,Y,'fmm',1,'verbose',0);

[usurfold,vsurfold] = evaluate_velocity_on_surface(solution_old,solution_old);
[usurfnew,vsurfnew] = evaluate_velocity_on_surface(solution_new,solution_new);

%% plot
figure('DefaultAxesFontSize',16);
contourf(X,Y,log10(abs(Uold-Unew)+eps),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
%caxis([-16,-12]);
axis equal;
title('$u$: $\log_{10}$(absolute error) old vs. new','interpreter','latex','fontsize',16);
%saveas(gca,'/home/david/Documents/Research/kernels/slip_tests/pipe/old_vs_new_matvec_u_err_amp_1e-1_eta_1_sq_matlab','epsc');

figure('DefaultAxesFontSize',16);
contourf(X,Y,log10(abs(Vold-Vnew)+eps),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
%caxis([-16,-12]);
axis equal;
title('$v$: $\log_{10}$(absolute error) old vs. new','interpreter','latex','fontsize',16);
%saveas(gca,'/home/david/Documents/Research/kernels/slip_tests/pipe/old_vs_new_matvec_v_err_amp_1e-1_eta_1_sq_matlab','epsc');

figure('DefaultAxesFontSize',16);
subplot(1,2,1);
semilogy(abs(usurfold));
xlim([1,N]);
grid on;
title('$u$ surf: abs err (old)','interpreter','latex','fontsize',16);

subplot(1,2,2);
semilogy(abs(vsurfold));
xlim([1,N]);
grid on;
title('$v$ surf: abs err (old)','interpreter','latex','fontsize',16);
%saveas(gca,'/home/david/Documents/Research/kernels/slip_tests/pipe/u_surf_err_old_amp_1e-1_eta_1_sq_matlab','epsc');

figure('DefaultAxesFontSize',16);
subplot(1,2,1);
semilogy(abs(usurfnew));
xlim([1,N]);
grid on;
title('$u$ surf: abs err (new)','interpreter','latex','fontsize',16);

subplot(1,2,2);
semilogy(abs(vsurfnew));
xlim([1,N]);
grid on;
title('$v$ surf: abs err (new)','interpreter','latex','fontsize',16);
%saveas(gca,'/home/david/Documents/Research/kernels/slip_tests/pipe/u_surf_err_new_amp_1e-1_eta_1_sq_matlab','epsc');
