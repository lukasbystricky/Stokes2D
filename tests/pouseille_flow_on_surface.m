% Pressure driven flow through a flat pipe, i.e. Pouseille flow. The exact
% solutions are known. Computes the on-surface values and compares them
% with the analytical ones for the following quantities:
%   - velocity
%   - pressure
%   - velocity gradient
%   - pressure gradient
%   - vorticity
%   - stress

%%%% Notes %%%%
% DK: 2022-07-07
% Accurate on-surface quantities for both a single-layer and combined
% formulation for no slip BC (alpha=0). When alpha=-1e-1 the single-layer
% formulation still yield very accurate results. However, the combined
% formulation has trouble converging and the error in some quantities are
% increased a couple of orders of magnitude. This holds true when using the
% velocity gradients to directly compute the slip BC.

close all;
clearvars;
clc;
format long;
set(groot,'defaultAxesTickLabelInterpreter','latex');

% create input structure
input_params = default_input_params('pouseuille_demo', 1);

% modify structure as needed
input_params.box_size = [1,2];
input_params.h = 0.5;    % pipe walls at +-0.5
input_params.panels = 5;
input_params.eta = inf;
input_params.plot_domain = 0;
input_params.alpha = -0.4;
input_params.slip = 1;
input_params.gmres_tol = 1e-12;
input_params.nbr_neighbor_pts = 4;

problem = flat_pipe_periodic(input_params);

z = problem.domain.z;
x = real(z);
y = imag(z);
N = length(z);

%% solve the problem
rhs = [zeros(N/2,1); zeros(N/2,1); zeros(N/2,1); zeros(N/2,1); problem.pressure_gradient_x; problem.pressure_gradient_y];
solution = solve_stokes(problem,rhs,'fmm',0);
solution.local_indices = 1:length(solution.q);

%% exact solutions
p = problem.pressure_gradient_x;
h = problem.h;
alpha = problem.alpha;

exact_solution_u = @(x,y) p/2*(y.^2-h^2) + h*p*alpha;
exact_solution_uy = @(x,y) p*y;
exact_solution_p = @(x,y) p*x;
exact_solution_dP = @(x,y) p;
exact_solution_omega = @(x,y) -p*y;

%% velocity
[Usurf, Vsurf] = evaluate_velocity_on_surface(solution, solution);
Usurf_exact = exact_solution_u(x,y);
Vsurf_exact = 0;

figure;
semilogy(abs(Usurf-Usurf_exact));
hold on;
semilogy(abs(Vsurf-Vsurf_exact));
grid on;
title('on-surface velocity','interpreter','latex','fontsize',16);
legend('$u$','$v$','interpreter','latex','fontsize',16);
ylabel('absolute error','interpreter','latex','fontsize',16);
xlabel('$n$','interpreter','latex','fontsize',16);
xlim([1,N]);

%% pressure
%solution.problem.domain.nbr_neighbor_pts = 4;
Psurf = evaluate_pressure_on_surface(solution, solution, 'fluid');
Psurf_exact = exact_solution_p(x,y);

figure;
semilogy(abs(Psurf-Psurf_exact)/max(abs(Psurf_exact))+eps);
grid on;
title('on-surface pressure','interpreter','latex','fontsize',16);
ylabel('relative error','interpreter','latex','fontsize',16);
xlabel('$n$','interpreter','latex','fontsize',16);
xlim([1,N]);

%% velocity gradient
[Uxsurf, Uysurf, Vxsurf, Vysurf] = evaluate_velocity_gradient_on_surface(solution, solution, 'fluid');
Uxsurf_exact = 0;
Uysurf_exact = exact_solution_uy(x,y);
Vxsurf_exact = 0;
Vysurf_exact = 0;

figure;
semilogy(abs(Uxsurf-Uxsurf_exact));
hold on;
semilogy(abs(Uysurf-Uysurf_exact)/max(abs(Uysurf_exact))+eps);
semilogy(abs(Vxsurf-Vxsurf_exact));
semilogy(abs(Vysurf-Vysurf_exact));
grid on;
title('on-surface velocity gradient','interpreter','latex','fontsize',16);
legend('$u_x$','$u_y$','$v_x$','$v_y$','interpreter','latex','fontsize',16);
ylabel('error','interpreter','latex','fontsize',16);
xlabel('$n$','interpreter','latex','fontsize',16);
xlim([1,N]);

%% pressure gradient
[Pxsurf, Pysurf] = evaluate_pressure_gradient_on_surface(solution, solution, 'fluid');
Pxsurf_exact = exact_solution_dP(x,y);
Pysurf_exact = 0;

figure;
semilogy(abs(Pxsurf-Pxsurf_exact)/max(abs(Pxsurf_exact))+eps);
hold on;
semilogy(abs(Pysurf-Pysurf_exact));
grid on;
title('on-surface pressure gradient','interpreter','latex','fontsize',16);
legend('$p_x$','$p_y$','interpreter','latex','fontsize',16);
ylabel('error','interpreter','latex','fontsize',16);
xlabel('$n$','interpreter','latex','fontsize',16);
xlim([1,N]);

%% vorticity
omegasurf = evaluate_vorticity_on_surface(solution, solution, 'fluid');
omegasurf_exact = exact_solution_omega(x,y);

figure;
semilogy(abs(omegasurf-omegasurf_exact)/max(abs(omegasurf_exact))+eps);
grid on;
title('on-surface vorticity','interpreter','latex','fontsize',16);
ylabel('relative error','interpreter','latex','fontsize',16);
xlabel('$n$','interpreter','latex','fontsize',16);
xlim([1,N]);

%% stress
%solution.problem.domain.nbr_neighbor_pts = 16;
[Sxxsurf, Sxysurf, Syxsurf, Syysurf] = evaluate_stress_on_surface(solution, solution, 'fluid');
Ssurf_exact = exact_solution_sigma_surface(problem);
Sxxsurf_exact = Ssurf_exact(:,1);
Syxsurf_exact = Ssurf_exact(:,2);
Sxysurf_exact = Ssurf_exact(:,3);
Syysurf_exact = Ssurf_exact(:,4);

figure;
semilogy(abs(Sxxsurf-Sxxsurf_exact)/max(abs(Sxxsurf_exact))+eps,'o-');
hold on;
semilogy(abs(Sxysurf-Sxysurf_exact)/max(abs(Sxysurf_exact))+eps,'x-');
semilogy(abs(Syxsurf-Syxsurf_exact)/max(abs(Syxsurf_exact))+eps,'d-');
semilogy(abs(Syysurf-Syysurf_exact)/max(abs(Syysurf_exact))+eps,'p-');
legend('$\sigma_{xx}$','$\sigma_{xy}$','$\sigma_{yx}$','$\sigma_{yy}$','interpreter','latex','fontsize',16);
grid on;
title('on-surface stress','interpreter','latex','fontsize',16);
ylabel('relative error','interpreter','latex','fontsize',16);
xlabel('$n$','interpreter','latex','fontsize',16);
xlim([1,N]);

%% more exact solutions
function S = exact_solution_sigma(X,Y,p)
[m,n] = size(X);
S = zeros(m,n,4);
for i = 1:m
    for j = 1:n
        x = X(i,j);
        y = Y(i,j);

        s = [-p*x p*y;
             p*y -p*x];

        S(i,j,:) = s(:);
    end
end
end

function sigma_surf = exact_solution_sigma_surface(problem)
p = problem.pressure_gradient_x;
z = problem.domain.z;

N = length(z);
sigma_surf = zeros(N,4);
for i = 1:N
    x = real(z(i));
    y = imag(z(i));
    s = exact_solution_sigma(x,y,p);
    sigma_surf(i,1) = s(1);
    sigma_surf(i,2) = s(2);
    sigma_surf(i,3) = s(3);
    sigma_surf(i,4) = s(4);
end
end
