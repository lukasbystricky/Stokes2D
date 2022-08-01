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
format long

% create input structure
input_params = default_input_params('pouseuille_demo', 1);

% modify structure as needed
input_params.box_size = [2,2];
input_params.h = 0.5;    % pipe walls at +-0.5
input_params.panels = 10;
input_params.eta = 1;
input_params.plot_domain = 0;
input_params.alpha = 0;
input_params.slip = 0;
%input_params.gmres_tol = 1e-8;

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
volume = 2*Lx*h;

exact_solution_u = @(x,y) p/2*(y.^2-h^2) + h*p*alpha;
exact_solution_v = @(x,y) zeros(size(x));
exact_solution_uy = @(x,y) p*y;
exact_solution_p = @(x,y) p*x;
exact_solution_dP = @(x,y) p;
exact_solution_omega = @(x,y) -p*y;

exact_solution_u_avg = 2*Lx*p*h^2*(alpha-h/3)/volume;
exact_solution_u_grad_avg = [0 0; 0 0];
exact_solution_p_avg = 0;
exact_solution_p_grad_avg = [(2*h*Lx*p)/volume; 0];
exact_solution_omega_avg = 0;
exact_solution_sigma_avg = [0 0; 0 0];

%% compute quantities
[Uc, Vc, X, Y, U, V] = evaluate_velocity(solution, 200, 'fmm', 0, 'verbose', 0);
[Pc, P, ~, ~] = evaluate_pressure(solution, X, Y, 'fmm', 0, 'verbose', 0);
[Uxc, Uyc, Vxc, Vyc, Ux, Uy, Vx, Vy] = evaluate_velocity_gradient(solution, X, Y);
[Px, Py] = evaluate_pressure_gradient(solution, X, Y);
[omegac, omega] = evaluate_vorticity(solution, X, Y);
[sxxc, sxyc, syxc, syyc, sxx, sxy, syx, syy, X, Y] = evaluate_stress(solution, X, Y);
sigma = exact_solution_sigma(X,Y,p);

% Stresses from pressure and velocity gradient calculations (works)
% velsxxc = -Pc+2*Uxc;
% velsxyc = Vxc+Uyc;
% velsyxc = Uyc+Vxc;
% velsyyc = -Pc+2*Vyc;

%% compute averages
solution.trim = 0;      % do not remove values outside the domain
[u_avg_alt1, u_avg_alt2, u_avg_alt3, u_grad_avg, p_avg, p_grad_avg, omega_avg, sigma_avg] = compute_pipe_averages(solution);
u_avg_alt1 = u_avg_alt1(1) + 1i*u_avg_alt1(2);
u_avg_alt2 = u_avg_alt2(1) + 1i*u_avg_alt2(2);
u_avg_alt3 = u_avg_alt3(1) + 1i*u_avg_alt3(2);

%% print averages errors
u_avg_alt1_err = abs(exact_solution_u_avg-u_avg_alt1)
u_avg_alt2_err = abs(exact_solution_u_avg-u_avg_alt2)
u_avg_alt3_err = abs(exact_solution_u_avg-u_avg_alt3)
u_grad_avg_err = abs(exact_solution_u_grad_avg-u_grad_avg)
p_avg_err = abs(exact_solution_p_avg-p_avg)
p_grad_avg_err = abs(exact_solution_p_grad_avg-p_grad_avg)
omega_avg_err = abs(exact_solution_omega_avg-omega_avg)
sigma_avg_err = abs(exact_solution_sigma_avg-sigma_avg)

%% plot velocities
figure;
subplot(2,2,1)
contourf(X,Y,Uc,20);
colorbar
axis equal
title('u');

subplot(2,2,2)
contourf(X,Y, log10(abs((Uc - exact_solution_u(X,Y))./...
    max(max(abs(exact_solution_u(X,Y)))))+eps),20);
colorbar
%caxis([-16,-1]);
axis equal
title('u: log_{10}(relative error)');

subplot(2,2,3)
contourf(X,Y,Vc,20);
colorbar
axis equal
title('v');

subplot(2,2,4)
contourf(X,Y, log10(abs(Vc - exact_solution_v(X,Y))+eps),20);
colorbar
%caxis([-16,-1]);
axis equal
title('v: log_{10}(relative error)');

%% plot velocity gradients
figure;
subplot(2,2,1)
contourf(X,Y,Uxc);
colorbar
axis equal
title('du/dx');

subplot(2,2,2)
contourf(X,Y,Uyc);
colorbar
axis equal
title('du/dy');

subplot(2,2,3)
contourf(X,Y,Vxc);
colorbar
axis equal
title('dv/dx');

subplot(2,2,4)
contourf(X,Y,Vyc);
colorbar
axis equal
title('dv/dy');

% errors
figure;
subplot(2,2,1)
contourf(X,Y, log10(abs(Uxc)+eps));
colorbar
%caxis([-16,-1]);
axis equal
title('dudx: log_{10}(absolute error)');

subplot(2,2,2)
contourf(X,Y, log10(abs((Uyc - exact_solution_uy(X,Y))./...
    max(max(abs(exact_solution_uy(X,Y)))))+eps));
colorbar
%caxis([-16,-1]);
axis equal
title('dudy: log_{10}(relative error)');

subplot(2,2,3)
contourf(X,Y, log10(abs(Vxc)+eps));
colorbar
%caxis([-16,-1]);
axis equal
title('dvdx: log_{10}(absolute error)');

subplot(2,2,4)
contourf(X,Y, log10(abs(Vxc)+eps));
colorbar
%caxis([-16,-1]);
axis equal
title('dudy: log_{10}(absolute error)');

%% plot pressure
figure;
subplot(1,2,1)
contourf(X,Y,Pc);
colorbar
axis equal
title('P');

subplot(1,2,2)
contourf(X,Y, log10(abs((Pc - exact_solution_p(X,Y))./...
    max(max(abs(exact_solution_uy(X,Y)))))+eps));
colorbar
%caxis([-16,-1]);
axis equal
title('P: log_{10}(relative error)');

%% plot pressure gradient
figure;
subplot(1,2,1)
contourf(X,Y,abs(Px + 1i*Py));
colorbar
axis equal
title('nabla p');

subplot(1,2,2)
contourf(X,Y,log10(abs(Px + 1i*Py  - exact_solution_dP(X,Y))+eps));
colorbar
%caxis([-16,-1]);
axis equal
title('nabla p: log_{10}(relative error)');

%% plot vorticity
subplot(1,2,1)
contourf(X,Y,omegac);
colorbar
axis equal
title('omega');

subplot(1,2,2)
contourf(X,Y, log10(abs((omegac - exact_solution_omega(X,Y))./...
    max(max(abs(exact_solution_omega(X,Y)))))+eps));
colorbar
%caxis([-16,-1]);
axis equal
title('omega: log_{10}(relative error)');

%% plot stress
figure;
subplot(2,2,1);
hold on;
contourf(X,Y,sxxc);
colorbar;
axis equal;
title('sxx');

subplot(2,2,2);
hold on;
contourf(X,Y,sxyc);
colorbar;
axis equal;
title('sxy');

subplot(2,2,3);
hold on;
contourf(X,Y,syxc);
colorbar;
axis equal;
title('syx');

subplot(2,2,4);
hold on;
contourf(X,Y,syyc);
colorbar;
axis equal;
title('syy');

% errors
figure;
subplot(2,2,1);
hold on;
contourf(X,Y,log10(abs(sxxc - sigma(:,:,1))./...
    max(max(abs(sigma(:,:,1)))) + eps));
colorbar;
axis equal;
title('sxx: log_{10} relative error');

subplot(2,2,2);
hold on;
contourf(X,Y,log10(abs(sxyc - sigma(:,:,3))./...
    max(max(abs(sigma(:,:,3)))) + eps));
colorbar;
axis equal;
title('sxy: log_{10} relative error');

subplot(2,2,3);
hold on;
contourf(X,Y,log10(abs(syxc - sigma(:,:,2))./...
    max(max(abs(sigma(:,:,2)))) + eps));
colorbar;
axis equal;
title('syx: log_{10} relative error');

subplot(2,2,4);
hold on;
contourf(X,Y,log10(abs(syyc - sigma(:,:,4))./...
    max(max(abs(sigma(:,:,4)))) + eps));
colorbar;
axis equal;
title('syy: log_{10} relative error');

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