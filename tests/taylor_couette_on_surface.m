% Taylor-Couette flow. Computes the on-surface values and compares them
% with the known analytical ones for the following quantities:
%   - velocity
%   - stress
% The slip length alpha must be set to zero for analytical solutions to be
% correct.

close all;
clearvars;
clc;
format long;
set(groot,'defaultAxesTickLabelInterpreter','latex');

% create input structure
input_params = default_input_params('couette', 0);

% modify structure as needed, or add additional problem-dependent params
input_params.panels = [10, 10];
input_params.plot_domain = 0;

% set up radii and centers, centers given as complex numbers (x+iy)
input_params.radii = [2, 0.5]; % the outer wall must be the first radii
input_params.centers = [0, 0];
input_params.omega = 1; % angular velocity of inner cylinder
input_params.alpha = 0;
input_params.slip = 0;

% initialize problem
problem = couette(input_params);
problem.domain.nbr_neighbor_pts = 16;

% domain
z = problem.domain.z;
x = real(z);
y = imag(z);
N = length(z);

% quantities for exact solution
ro = input_params.radii(1);
ri = input_params.radii(2);
omega = input_params.omega;
alpha = -input_params.alpha;    % note minus sign
Aold = -omega/(ro^2*(1/ri^2-1/ro^2));
Bold = omega/(1/ri^2-1/ro^2);
Anew = (1-2*alpha/ro)/ro^2;
Bnew = (1+2*alpha/ri)/ri^2;

% exact solution in polar coordinates
exact_solution_r = @(x,y) zeros(size(x));
exact_solution_theta = @(x,y) (omega/(Anew-Bnew)) * (Anew*sqrt(x.^2 + y.^2) - 1./sqrt(x.^2 + y.^2));
sigma_exact = exact_solution_sigma(x,y,Aold,Bold);

%% solve the problem
rhs = [zeros(N/2,1);
       ri*omega*ones(N/2,1);
       zeros(N,1)];
solution = solve_stokes(problem,rhs,'fmm',0);
solution.local_indices = 1:length(solution.q);

%% velocity
[u,v] = evaluate_velocity_on_surface(solution,solution);
ur = (u.*x + v.*y)./sqrt(x.^2 + y.^2);
utheta = (-u.*y + v.*x)./sqrt(x.^2 + y.^2);

ur_exact = exact_solution_r(x,y);
utheta_exact = exact_solution_theta(x,y);

figure;
semilogy(abs(ur-ur_exact));
hold on;
semilogy(abs(utheta-utheta_exact));
grid on;
title('on-surface velocity','interpreter','latex','fontsize',16);
legend('$u_r$','$u_\theta$','interpreter','latex','fontsize',16);
ylabel('absolute error','interpreter','latex','fontsize',16);
xlabel('$n$','interpreter','latex','fontsize',16);
xlim([1,N]);

%% stress
[sxx,sxy,syx,syy] = evaluate_stress_on_surface(solution,solution,'fluid');

sxx_exact = sigma_exact(:,1);
syx_exact = sigma_exact(:,2);
sxy_exact = sigma_exact(:,3);
syy_exact = sigma_exact(:,4);

figure;
semilogy(abs(sxx-sxx_exact)/max(abs(sxx_exact))+eps,'p-');
hold on;
semilogy(abs(sxy-sxy_exact)/max(abs(sxy_exact))+eps,'x-');
semilogy(abs(syx-syx_exact)/max(abs(syx_exact))+eps,'d-');
semilogy(abs(syy-syy_exact)/max(abs(syy_exact))+eps,'p-');
legend('$\sigma_{xx}$','$\sigma_{xy}$','$\sigma_{yx}$','$\sigma_{yy}$','interpreter','latex','fontsize',16);
grid on;
title('on-surface stress','interpreter','latex','fontsize',16);
ylabel('relative error','interpreter','latex','fontsize',16);
xlabel('$n$','interpreter','latex','fontsize',16);
xlim([1,N]);

%% functions
function S = exact_solution_sigma(X,Y,A,B)
[m,n] = size(X);
S = zeros(m,n,4);
for i = 1:m
    for j = 1:n
        x = X(i,j);
        y = Y(i,j);

        exact_solution_ux = (2*B*x.*y)./(x.^2+y.^2).^2;
        exact_solution_uy = B*(y.^2-x.^2)./(x.^2+y.^2).^2 - A;
        exact_solution_vx = B*(y.^2-x.^2)./(x.^2+y.^2).^2 + A;
        exact_solution_vy = -(2*B*x.*y)./(x.^2+y.^2).^2;

        U = [exact_solution_ux exact_solution_vx;
             exact_solution_uy exact_solution_vy];
        s = U + U';
        S(i,j,:) = s(:);
    end
end
end