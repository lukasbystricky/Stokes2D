% Check consistency of the solution by adding two reference cells in the
% periodic direction.

%%%% Notes %%%%
% DK: 2022-07-07
% For alpha=-1e-1 and amp=1e-1, then the test works only if we use the
% velocity gradients to compute the slip BC.

close all;
clearvars;
clc;
set(groot,'defaultAxesTickLabelInterpreter','latex');

% parameters
alpha = -1e-1;
amp = 1e-1;

% create input structure
input_params = default_input_params('slip_pipe_test', 1);

% modify structure as needed
input_params.box_size = [1,2];
input_params.h = 0.5;
input_params.panels = 10;
input_params.plot_domain = 0;
input_params.slip = 1;
input_params.alpha = alpha;
input_params.A = amp;
input_params.d = 0.5;
input_params.eta = inf;
input_params.gmres_tol = 1e-10;

% smaller problem
problem1 = sine_pipe_periodic(input_params);

% problem with added periodicity at both left and right side. Note that
% since the reference cell is three times as large in the x-direction, we
% must compensate the pressure gradient in order for the solutions to be
% comparable
input_params.box_size = [3,2];
input_params.pressure_drop_x = 3*input_params.pressure_drop_x;
input_params.panels = 30;
problem2 = sine_pipe_periodic(input_params);

% plot domains
figure('DefaultAxesFontSize',16);
hold on;
plot(problem2.domain.z,'x')
plot(problem1.domain.z,'o')
xlabel('x');
ylabel('y');
xlim([-3/2,3/2]);
grid on;
legend('problem 2','problem 1');

%% solve problems
N1 = length(problem1.domain.z);
rhs1 = [zeros(2*N1,1); problem1.pressure_gradient_x; problem1.pressure_gradient_y];
solution1 = solve_stokes(problem1,rhs1,'fmm',0);
solution1.local_indices = 1:length(solution1.q);

N2 = length(problem2.domain.z);
rhs2 = [zeros(2*N2,1); problem2.pressure_gradient_x; problem2.pressure_gradient_y];
solution2 = solve_stokes(problem2,rhs2,'fmm',0);
solution2.local_indices = 1:length(solution2.q);

%% compute velocities
[U1,V1,X,Y] = evaluate_velocity(solution1,200,'fmm',0,'verbose',0);
[Usurf1,Vsurf1] = evaluate_velocity_on_surface(solution1,solution1);
P1 = evaluate_pressure(solution1,X,Y,'fmm',0,'verbose',0);

[U2,V2] = evaluate_velocity(solution2,X,Y,'fmm',0,'verbose',0);
[Usurf2,Vsurf2] = evaluate_velocity_on_surface(solution2,solution2);
P2 = evaluate_pressure(solution2,X,Y,'fmm',0,'verbose',0);

%% plot velocities
figure('DefaultAxesFontSize',16);
subplot(2,2,1);
contourf(X,Y,U1,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$u_1$','interpreter','latex','fontsize',16);

subplot(2,2,2);
contourf(X,Y,V1,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$v_1$','interpreter','latex','fontsize',16);

subplot(2,2,3);
contourf(X,Y,U2,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$u_2$','interpreter','latex','fontsize',16);

subplot(2,2,4);
contourf(X,Y,V2,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$v_2$','interpreter','latex','fontsize',16);

%% plot velocity errors
figure('DefaultAxesFontSize',16);
subplot(1,2,1);
contourf(X,Y,log10(abs(U1-U2)+eps),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$u$ $\log_{10}$(absolute error)','interpreter','latex','fontsize',16);

subplot(1,2,2);
contourf(X,Y,log10(abs(V1-V2)+eps),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$v$ $\log_{10}$(absolute error)','interpreter','latex','fontsize',16);

%% plot pressures and errors
figure('DefaultAxesFontSize',16);
subplot(2,2,1);
contourf(X,Y,P1,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$p_1$','interpreter','latex','fontsize',16);

subplot(2,2,2);
contourf(X,Y,P2,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('$p_2$','interpreter','latex','fontsize',16);

subplot(2,2,3);
contourf(X,Y,log10(abs(P1-P2)+eps),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
title('p absolute error','interpreter','latex','fontsize',16);

%% plot on-surface velocity
figure('DefaultAxesFontSize',16);
subplot(2,2,1);
plot(Usurf1);
xlim([1,length(Usurf1)]);
grid on;
title('on-surface $u_1$','interpreter','latex','fontsize',16);

subplot(2,2,2);
plot(Vsurf1);
xlim([1,length(Vsurf1)]);
grid on;
title('on-surface $v_1$','interpreter','latex','fontsize',16);

subplot(2,2,3);
plot(Usurf2);
xlim([1,length(Usurf2)]);
grid on;
title('on-surface $u_2$','interpreter','latex','fontsize',16);

subplot(2,2,4);
plot(Vsurf2);
xlim([1,length(Vsurf2)]);
grid on;
title('on-surface $v_2$','interpreter','latex','fontsize',16);

%% plot density error
idx_top = N2/6+1:N2/6+N1/2;
idx_bottom = 2*N2/3+1:2*N2/3+N1/2;
qerr_real = abs(solution1.q(:,1)-[solution2.q(idx_top,1); solution2.q(idx_bottom,1)]);
qerr_imag = abs(solution1.q(:,2)-[solution2.q(idx_top,2); solution2.q(idx_bottom,2)]);

figure;
subplot(1,2,1);
semilogy(qerr_real);
title('$Re(q)$','interpreter','latex','fontsize',16);
xlabel('discretization points');
ylabel('absolute error');
grid on;

subplot(1,2,2);
semilogy(qerr_imag);
title('$Im(q)$','interpreter','latex','fontsize',16);
xlabel('discretization points');
ylabel('absolute error');
grid on;
