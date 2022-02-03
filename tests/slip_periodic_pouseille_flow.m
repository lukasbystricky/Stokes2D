% Code to test the slip boundary condition on flat and curved boundaries by
% doing the following:
%
% 1. Let alpha=0 and solve pressure driven flow with non-zero boundary
%    condition.
% 2. Use solution of (1) and read of new boundary condition, e.g. finding
%    new alpha values.
% 3. Set the boundary conditions of (2) and solve the new problem using
%    both GMRES and by performing a direct solve (DS).
% 4. The solution from (1) and (3) should now be the same.

close all;
clearvars;
clc;

% create input structure
input_params = default_input_params('pouseuille_demo', 1);

% modify structure as needed
input_params.box_size = [4,6];
input_params.h = 0.5;    % pipe walls at +-0.5
input_params.panels = 10;
input_params.eta = 1;
input_params.plot_domain = 0;
input_params.alpha = 0;
input_params.A = 5*1e-2;

%problem = flat_pipe_periodic(input_params);
problem = sine_pipe_periodic(input_params);
Lx = problem.Lx;
Ly = problem.Ly;

% solve the problem with no slip
N = length(problem.domain.z);
rhs_ns = [0.001*ones(N/2,1); 0.02*ones(N/2,1); 0.03*zeros(N/2,1); 0.01*zeros(N/2,1); problem.pressure_gradient_x; problem.pressure_gradient_y];
%rhs_ns = [zeros(2*N,1); problem.pressure_gradient_x; problem.pressure_gradient_y];
solution_ns = solve_stokes(problem,rhs_ns,'fmm',0);
solution_ns.local_indices = 1:length(solution_ns.q);

%% compute alpha
domain = problem.domain;

% normalized tangent
nu1 = real(domain.zp)./abs(domain.zp);
nu2 = imag(domain.zp)./abs(domain.zp);

% normalized normal
n1 = real(-1i*domain.zp)./abs(domain.zp);
n2 = imag(-1i*domain.zp)./abs(domain.zp);

% on-surface values
solution_ns.local_indices = 1:length(solution_ns.q);

[u1,u2] = evaluate_velocity_on_surface(solution_ns, solution_ns);
udotn = u1.*n1 + u2.*n2;
udotnu = u1.*nu1 + u2.*nu2;

[sxx,sxy,syx,syy] = evaluate_stress_on_surface(solution_ns, solution_ns, 'fluid');
t1 = n1.*sxx + n2.*syx;
t2 = n1.*sxy + n2.*syy;
nudott = t1.*nu1 + t2.*nu2;

% [ux,uy,vx,vy] = evaluate_velocity_gradient_on_surface(solution_ns, solution_ns, 'fluid');
% nudott = (vx+uy).*n2.*nu1;

alpha = -udotnu./nudott;
problem.alpha = alpha;
max(abs(alpha))

%% solve with slip length alpha (GMRES)
rhs_s = [zeros(N/2,1); zeros(N/2,1); udotn; problem.pressure_gradient_x; problem.pressure_gradient_y];
solution_s_gmres = solve_stokes(problem,rhs_s,'fmm',0);
solution_s_gmres.local_indices = 1:length(solution_s_gmres.q);

%% solve with slip length alpha (direct)
rhs_s = [zeros(N/2,1); zeros(N/2,1); udotn; problem.pressure_gradient_x; problem.pressure_gradient_y];
[solution_s_ds,K_s] = solve_direct(problem,rhs_s);
size(K_s)
rank(K_s)
cond(K_s)

%% evaluate velocity
[U_ns, V_ns, X, Y, ~, ~] = evaluate_velocity(solution_ns, 200, 'fmm', 1, 'verbose', 0);
[U_s_gmres, V_s_gmres, ~, ~, ~, ~] = evaluate_velocity(solution_s_gmres, 200, 'fmm', 1, 'verbose', 0);
[U_s_ds, V_s_ds, ~, ~, ~, ~] = evaluate_velocity(solution_s_gmres, 200, 'fmm', 1, 'verbose', 0);

[Usurf_ns, Vsurf_ns] = evaluate_velocity_on_surface(solution_ns, solution_ns);
[Usurf_s_gmres, Vsurf_s_gmres] = evaluate_velocity_on_surface(solution_s_gmres, solution_s_gmres);
[Usurf_s_ds, Vsurf_s_ds] = evaluate_velocity_on_surface(solution_s_ds, solution_s_ds);

%% plot
plot_gmres = 1;

if plot_gmres
    U_s = U_s_gmres;
    V_s = V_s_gmres;
    Usurf_s = Usurf_s_gmres;
    Vsurf_s = Vsurf_s_gmres;
    solution_s = solution_s_gmres;
else
    U_s = U_s_ds;
    V_s = V_s_ds;
    Usurf_s = Usurf_s_ds;
    Vsurf_s = Vsurf_s_ds;
    solution_s = solution_s_ds;
end

figure;
subplot(3,2,1)
contourf(X,Y,U_ns,20);
colorbar
axis equal
title('u no-slip');

subplot(3,2,2)
contourf(X,Y,V_ns,20);
colorbar
axis equal
title('v no-slip');

subplot(3,2,3)
contourf(X,Y,U_s,20);
colorbar
axis equal
title('u slip');

subplot(3,2,4)
contourf(X,Y,V_s,20);
colorbar
axis equal
title('v slip');

subplot(3,2,5)
contourf(X,Y,log10(abs(U_ns-U_s)./max(max(abs(U_ns)))+eps),20);
colorbar
axis equal
title('u: log_{10}(relative error)');

subplot(3,2,6)
contourf(X,Y,log10(abs(V_ns-V_s)+eps),20);
colorbar
axis equal
title('v: log_{10}(absolute error)');

figure;
semilogy(abs(Usurf_ns-Usurf_s)+eps);
hold on;
semilogy(abs(Vsurf_ns-Vsurf_s)+eps);
grid on;
legend('u','v');
title('on-surface velocity');
xlabel('discretization points');
ylabel('absolute error');

figure;
semilogy(abs(solution_ns.q-solution_s.q));
grid on;
norm(solution_ns.q(:,1)-solution_s.q(:,1))
norm(solution_ns.q(:,2)-solution_s.q(:,2))
legend('real','imag');
title('slip vs no-slip density');

disp('mean velocities');
solution_ns.u_avg
solution_s.u_avg

function [solution,K] = solve_direct(problem,rhs)
N = length(problem.domain.z);
K = zeros(2*N+2,2*N+2);
XX = eye(2*N+2);

parfor i = 1:2*N+2
    K(:,i) = matvec_combined_robin(XX(:,i),problem);
end

X = K\rhs;
solution.q = [X(1:N), X(N+1:2*N)];
solution.u_avg = [X(end-1); X(end)];
solution.trim = 1;
solution.local_indices = 1:length(solution.q);
solution.problem = problem;
end
