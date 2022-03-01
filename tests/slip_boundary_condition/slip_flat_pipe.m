clearvars
clc
close all
format long
set(groot,'defaultAxesTickLabelInterpreter','latex');

% create input structure
input_params = default_input_params('pouseuille_demo', 1);

% modify structure as needed
input_params.box_size = [2,5];
input_params.h = 0.5;    % pipe walls at +-0.5
input_params.panels = 10;
input_params.eta = inf;
input_params.plot_domain = 0;
input_params.alpha = 0;
input_params.slip = 1;
input_params.d = 0.5;

problem = flat_pipe_periodic(input_params);
Lx = problem.Lx;
Ly = problem.Ly;

%% solve the problem
N = length(problem.domain.z);
rhs = [zeros(N/2,1); zeros(N/2,1); zeros(N/2,1); zeros(N/2,1); 
       problem.pressure_gradient_x; problem.pressure_gradient_y];
solution = solve_stokes(problem,rhs,'fmm',0);
solution.local_indices = 1:length(solution.q);

%% exact solution
p = problem.pressure_gradient_x;
h = problem.h;
alpha = problem.alpha;
x = real(problem.domain.z);
y = imag(problem.domain.z);

exact_solution_u = @(x,y) p/2*(y.^2-h^2) + h*p*alpha;
exact_solution_v = @(x,y) zeros(size(x));
sigma = exact_solution_sigma(x,y,p);

%% compute quantities
[Uc, Vc, X, Y, U, V] = evaluate_velocity(solution, 200, 'fmm', 0, 'verbose', 0);
[usurf,vsurf] = evaluate_velocity_on_surface(solution,solution);
[sxxsurf,sxysurf,syxsurf,syysurf] = evaluate_stress_on_surface(solution,solution,'fluid');

%% plot
figure('DefaultAxesFontSize',16);
subplot(1,2,1)
contourf(X,Y, log10(abs((Uc - exact_solution_u(X,Y))./...
    max(max(abs(exact_solution_u(X,Y)))))+eps),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
%caxis([-16,-1]);
axis equal
title('$u$: $\log_{10}$(relative error)','interpreter','latex','fontsize',16);

subplot(1,2,2)
contourf(X,Y, log10(abs(Vc - exact_solution_v(X,Y))+eps),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
%caxis([-16,-1]);
axis equal
title('$v$: $\log_{10}$(absolute error)','interpreter','latex','fontsize',16);
%saveas(gca,'/home/david/Documents/Research/kernels/slip_tests/pipe/flat_u_err_alpha_1e-1_eta_inf_sq_matlab','epsc');

%%
figure('DefaultAxesFontSize',16);
subplot(2,2,1);
semilogy(abs((sxxsurf-sigma(:,1,1))./max(sigma(:,1,1)))+eps);
title('$\sigma_{xx}$: relative error','interpreter','latex','fontsize',16);
grid on;
xlim([1,length(x)]);

subplot(2,2,2);
semilogy(abs((sxysurf-sigma(:,1,3))./max(sigma(:,1,3)))+eps);
title('$\sigma_{xy}$: relative error','interpreter','latex','fontsize',16);
grid on;
xlim([1,length(x)]);

subplot(2,2,3);
semilogy(abs((syxsurf-sigma(:,1,2))./max(sigma(:,1,2)))+eps);
title('$\sigma_{yx}$: relative error','interpreter','latex','fontsize',16);
grid on;
xlim([1,length(x)]);

subplot(2,2,4);
semilogy(abs((syysurf-sigma(:,1,4))./max(sigma(:,1,4)))+eps);
title('$\sigma_{yy}$: relative error','interpreter','latex','fontsize',16);
grid on;
xlim([1,length(x)]);
%saveas(gca,'/home/david/Documents/Research/kernels/slip_tests/pipe/flat_sigma_surf_err_alpha_1e-1_eta_inf_sq_matlab','epsc');

%% functions
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
