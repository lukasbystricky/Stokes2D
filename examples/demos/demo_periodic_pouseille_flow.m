% Basic demo script. Computes the velocity field for a pressure driven flow
% through a flat pipe, i.e. Pouseille flow. The exact solution is known and
% given by -y*(y-1)*dP/Lx/2, where dP/Lx is the imposed pressure gradient.
% Note that this problem is only periodic in the x direction, however we
% can still use the 2P periodic spectral Ewald code be embedding the 
% channel in a periodic box. We use a combined-layer formulation where the 
% pressure gradient is balanced by the single-layer potential. 

close all
clearvars
clc

% create input structure
input_params = default_input_params('pouseuille_demo', 1);

% modify structure as needed
input_params.box_size = [3,5];
input_params.h = 0.5;    % pipe walls at +-0.5
input_params.panels = 20;
input_params.eta = 1;
input_params.plot_domain = 0;

problem = flat_pipe_periodic(input_params);
Lx = problem.Lx;
Ly = problem.Ly;

%% solve the problem
solution = solve_stokes(problem,'fmm',0);

%% exact solution
p = problem.pressure_gradient_x;
h = problem.h;
volume = 2*Lx*h;

exact_solution_u = @(x,y) p/2*(y.^2-h^2);
exact_solution_uy = @(x,y) p*y;
exact_solution_p = @(x,y) p*x;
exact_solution_dP = @(x,y) p;
exact_solution_omega = @(x,y) -p*y;

exact_solution_u_avg = (-2*h^3*Lx*p)/3/volume;
exact_solution_u_grad_avg = [0 0; 0 0];
exact_solution_p_avg = 0;
exact_solution_p_grad_avg = [(2*h*Lx*p)/volume; 0];

%% compute quantities
[Uc, Vc, X, Y, U, V] = evaluate_velocity(solution, 100, 'fmm', 0, 'verbose', 0);
[Pc, P, ~, ~] = evaluate_pressure(solution, X, Y, 'fmm', 0, 'verbose', 0);
[Uxc, Uyc, Vxc, Vyc, Ux, Uy, Vx, Vy] = evaluate_velocity_gradient(solution, X, Y);
[Px, Py] = evaluate_pressure_gradient(solution, X, Y);
[omegac, omega] = evaluate_vorticity(solution, X, Y);

%% compute averages
solution.trim = 0;      % do not remove values outside the domain
[u_avg, u_grad_avg, p_avg, p_grad_avg] = compute_pipe_averages(solution);
u_avg = u_avg(1) + 1i*u_avg(2);

%% print averages errors
u_avg_err = abs(exact_solution_u_avg-u_avg)
u_grad_avg_err = abs(exact_solution_u_grad_avg-u_grad_avg)
p_avg_err = abs(exact_solution_p_avg-p_avg)
p_grad_avg_err = abs(exact_solution_p_grad_avg-p_grad_avg)

%% plot
figure;
subplot(5,2,1)
contourf(X,Y,Uc);
colorbar
axis equal
title('u');

subplot(5,2,2)
contourf(X,Y, log10(abs((Uc - exact_solution_u(X,Y))./...
    max(max(abs(exact_solution_u(X,Y)))))+eps));
colorbar
%caxis([-16,-1]);
axis equal
title('u: log_{10}(relative error)');

subplot(5,2,3)
contourf(X,Y,Uyc);
colorbar
axis equal
title('du/dy');

subplot(5,2,4)
contourf(X,Y, log10(abs((Uyc - exact_solution_uy(X,Y))./...
    max(max(abs(exact_solution_uy(X,Y)))))+eps));
colorbar
%caxis([-16,-1]);
axis equal
title('dudy: log_{10}(relative error)');

subplot(5,2,5)
contourf(X,Y,Pc);
colorbar
axis equal
title('P');

subplot(5,2,6)
contourf(X,Y, log10(abs((Pc - exact_solution_p(X,Y))./...
    max(max(abs(exact_solution_uy(X,Y)))))+eps));
colorbar
%caxis([-16,-1]);
axis equal
title('P: log_{10}(relative error)');

subplot(5,2,7)
contourf(X,Y,abs(Px + 1i*Py));
colorbar
axis equal
title('nabla p');

subplot(5,2,8)
contourf(X,Y,log10(abs(Px + 1i*Py  - exact_solution_dP(X,Y))+eps));
colorbar
%caxis([-16,-1]);
axis equal
title('nabla p: log_{10}(relative error)');

subplot(5,2,9)
contourf(X,Y,omegac);
colorbar
axis equal
title('omega');

subplot(5,2,10)
contourf(X,Y, log10(abs((omegac - exact_solution_omega(X,Y))./...
    max(max(abs(exact_solution_omega(X,Y)))))+eps));
colorbar
%caxis([-16,-1]);
axis equal
title('omega: log_{10}(relative error)');

%% individual figures
set(groot,'defaultAxesTickLabelInterpreter','latex');

figure;
contourf(X,Y, log10(abs((Uc - exact_solution_u(X,Y))./...
    max(max(abs(exact_solution_u(X,Y)))))+eps));
colorbar
caxis([-16,-10]);
axis equal;
xlabel({'$x$'},'interpreter','latex','fontsize',16);
ylabel({'$y$'},'interpreter','latex','fontsize',16);
title('Relative error $\mathbf{u}$ ($\log_{10}$)','interpreter','latex','fontsize',16);
%saveas(gca,'/afs/kth.se/home/d/a/davkra/Documents/phd_project_presentation_2021/rel_err_u','epsc');

figure;
contourf(X,Y, log10(abs((Pc - exact_solution_p(X,Y))./...
    max(max(abs(exact_solution_uy(X,Y)))))+eps));
colorbar
caxis([-16,-10]);
axis equal;
xlabel({'$x$'},'interpreter','latex','fontsize',16);
ylabel({'$y$'},'interpreter','latex','fontsize',16);
title('Relative error $p$ ($\log_{10}$)','interpreter','latex','fontsize',16);
%saveas(gca,'/afs/kth.se/home/d/a/davkra/Documents/phd_project_presentation_2021/rel_err_p','epsc');

figure;
contourf(X,Y, log10(abs((Uyc - exact_solution_uy(X,Y))./...
    max(max(abs(exact_solution_uy(X,Y)))))+eps));
colorbar
caxis([-16,-10]);
axis equal;
xlabel({'$x$'},'interpreter','latex','fontsize',16);
ylabel({'$y$'},'interpreter','latex','fontsize',16);
title('Relative error $\nabla\mathbf{u}$ ($\log_{10}$)','interpreter','latex','fontsize',16);
%saveas(gca,'/afs/kth.se/home/d/a/davkra/Documents/phd_project_presentation_2021/rel_err_nabla_u','epsc');

figure;
contourf(X,Y,log10(abs(Px + 1i*Py  - exact_solution_dP(X,Y))+eps));
colorbar
caxis([-16,-10]);
axis equal;
xlabel({'$x$'},'interpreter','latex','fontsize',16);
ylabel({'$y$'},'interpreter','latex','fontsize',16);
title('Absolute error $\nabla p$ ($\log_{10}$)','interpreter','latex','fontsize',16);
%saveas(gca,'/afs/kth.se/home/d/a/davkra/Documents/phd_project_presentation_2021/abs_err_nabla_p','epsc');
