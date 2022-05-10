% Test Dirichlet boundary conditions using the slip matvec code.
%   1. Solve pressure driven with flat walls with alpha=0.
%   2. Put a curved boundary inside the two flat walls.
%   3. Read of boundary condition on curved wall.
%   4. Solve smaller problem with curved wall (and alpha=0).
%   5. Compare solutions inside the domain.

close all;
clearvars;
clc;
set(groot,'defaultAxesTickLabelInterpreter','latex');

% parameters
alpha_flat = 1e-1;
alpha_curve = 0;
amp = 1e-1;
use_exact_bc = 1;

% create input structure
input_params = default_input_params('curved_pipe_slip_noslip_test', 1);

% modify structure as needed
input_params.box_size = [2,2];
input_params.h = 0.5;
input_params.panels = 10;
input_params.eta = inf;
input_params.plot_domain = 0;
input_params.slip = 0;

input_params.alpha = alpha_flat;
input_params.A = amp;
input_params.d = 0.3;

%problem_flat = flat_pipe_periodic(input_params);
input_params.d = 0.5;
input_params.A = 0;
problem_flat = sine_pipe_periodic(input_params);
input_params.A = amp;
input_params.d = 0.3;
input_params.slip = 1;
input_params.pressure_drop_x = 1;
problem_curve = sine_pipe_periodic(input_params);

% exact solutions
p = problem_flat.pressure_gradient_x;
h = problem_flat.h;
exact_solution_u = @(x,y) p/2*(y.^2-h^2) + h*p*alpha_flat;
exact_solution_v = @(x,y) zeros(size(x));

%% solve with flat walls
N = length(problem_flat.domain.z);
rhs_flat = [zeros(2*N,1); problem_flat.pressure_gradient_x; problem_flat.pressure_gradient_y];
solution_flat = solve_stokes(problem_flat,rhs_flat,'fmm',0);

%% read of boundary condition on upper curved and lower flat wall
domain_curve = problem_curve.domain;

% normalized tangent
nu1 = real(domain_curve.zp)./abs(domain_curve.zp);
nu2 = imag(domain_curve.zp)./abs(domain_curve.zp);

% normalized normal
n1 = real(-1i*domain_curve.zp)./abs(domain_curve.zp);
n2 = imag(-1i*domain_curve.zp)./abs(domain_curve.zp);

% points on curve
zcurve = problem_curve.domain.z(1:N/2);
xcurve = real(zcurve);
ycurve = imag(zcurve);

% points flat walls
zflat = problem_flat.domain.z;
xflat = real(zflat);
yflat = imag(zflat);

% get boundary values
if use_exact_bc
    % curve velocity
    u1curve = exact_solution_u(xcurve,ycurve);
    u2curve = zeros(N/2,1);
    
    % curve stress
    sigmacurve = exact_solution_sigma(xcurve,ycurve,p);
    sxxcurve = sigmacurve(:,1,1);
    syxcurve = sigmacurve(:,1,2);
    sxycurve = sigmacurve(:,1,3);
    syycurve = sigmacurve(:,1,4);
    
    % flat velocity
    u1flat = exact_solution_u(xflat,yflat);
    u2flat = zeros(N,1);
    
    % flat stress
    sigmaflat = exact_solution_sigma(xflat,yflat,p);
    sxxflat = sigmaflat(:,1,1);
    syxflat = sigmaflat(:,1,2);
    sxyflat = sigmaflat(:,1,3);
    syyflat = sigmaflat(:,1,4);
else
    % curve velocity and stress
    [u1curve,u2curve] = evaluate_velocity(solution_flat,xcurve,ycurve);
    [sxxcurve,sxycurve,syxcurve,syycurve] = evaluate_stress(solution_flat,xcurve,ycurve);

    % flat velocity and stress
    solution_flat.local_indices = 1:N;
    [u1flat,u2flat] = evaluate_velocity_on_surface(solution_flat, solution_flat);
    [sxxflat,sxyflat,syxflat,syyflat] = evaluate_stress_on_surface(solution_flat, solution_flat, 'fluid');
end

% combine boundaries
u1 = [u1curve; u1flat(N/2+1:N)];
u2 = [u2curve; u2flat(N/2+1:N)];
sxx = [sxxcurve; sxxflat(N/2+1:N)];
sxy = [sxycurve; sxyflat(N/2+1:N)];
syx = [syxcurve; syxflat(N/2+1:N)];
syy = [syycurve; syyflat(N/2+1:N)];

% normal direction
udotn = u1.*n1 + u2.*n2;

% tangential direction
udotnu = u1.*nu1 + u2.*nu2;
t1 = n1.*sxx + n2.*syx;
t2 = n1.*sxy + n2.*syy;
nudott = t1.*nu1 + t2.*nu2;

% boundary condition
problem_curve.alpha = [alpha_curve*ones(N/2,1); alpha_flat*ones(N/2,1)];
g1_curve = udotn(1:N/2);
g2_curve = udotnu(1:N/2) + problem_curve.alpha(1:N/2).*nudott(1:N/2);
g1_flat = udotn(N/2+1:end);
g2_flat = udotnu(N/2+1:end) + problem_curve.alpha(N/2+1:end).*nudott(N/2+1:end);
g1 = [g1_curve; g1_flat];
g2 = [g2_curve; g2_flat];

rhs_curve = [g2; g1; problem_curve.pressure_gradient_x; problem_curve.pressure_gradient_y];
%rhs_curve = [g2; g1; 0; 0];

%% solve new problem with curved top wall
solution_curve = solve_stokes(problem_curve,rhs_curve,'fmm',0);
solution_curve.local_indices = 1:length(solution_curve.q);
%[solution_curve,K] = solve_direct(problem_curve,rhs_curve);

%% compute velocities off-surface
[Ucurve,Vcurve,Xcurve,Ycurve] = evaluate_velocity(solution_curve,200,'fmm',0,'verbose',0);
[Uflat,Vflat,Xflat,Yflat] = evaluate_velocity(solution_flat,200,'fmm',0,'verbose',0);
[Uflat_inside,Vflat_inside] = evaluate_velocity(solution_flat,Xcurve,Ycurve,'fmm',0,'verbose',0);
Uexact = exact_solution_u(Xcurve,Ycurve);
Vexact = exact_solution_v(Xcurve,Ycurve);

ymin = min(Yflat(:));
ymax = max(Yflat(:));

%% compute velocities on-surface
[usurf,vsurf] = evaluate_velocity_on_surface(solution_curve,solution_curve);
[sxxsurf,sxysurf,syxsurf,syysurf] = evaluate_stress_on_surface(solution_curve,solution_curve,'fluid');

udotnsurf = usurf.*n1 + vsurf.*n2;
udotnusurf = usurf.*nu1 + vsurf.*nu2;
t1surf = n1.*sxxsurf + n2.*syxsurf;
t2surf = n1.*sxysurf + n2.*syysurf;
nudottsurf = t1surf.*nu1 + t2surf.*nu2;

g1surf = udotnsurf;
g2surf = udotnusurf + problem_curve.alpha.*nudottsurf;

%% plot solutions
figure('DefaultAxesFontSize',16);
subplot(2,2,1);
contourf(Xflat,Yflat,Uflat,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
caxis('auto');
lim = caxis;
axis equal;
ylim([ymin,ymax]);
title('$u$ flat','interpreter','latex','fontsize',16);

subplot(2,2,2);
contourf(Xflat,Yflat,Vflat,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
ylim([ymin,ymax]);
title('$v$ flat','interpreter','latex','fontsize',16);

subplot(2,2,3);
contourf(Xcurve,Ycurve,Ucurve,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
caxis(lim);
axis equal;
ylim([ymin,ymax]);
title('$u$ curve','interpreter','latex','fontsize',16);

subplot(2,2,4);
contourf(Xcurve,Ycurve,Vcurve,20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
axis equal;
ylim([ymin,ymax]);
title('$v$ curve','interpreter','latex','fontsize',16);
%saveas(gca,'/home/david/Documents/Research/kernels/slip_tests/pipe/vel_amp_1e-6_alpha_1e-5_eta_inf_sq_matlab','epsc');

%% plot errors
figure('DefaultAxesFontSize',16);
subplot(2,2,1);
contourf(Xcurve,Ycurve,log10(abs((Uexact-Uflat_inside)./max(max(abs(Uexact))))+eps),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
%caxis([-16,-12]);
axis equal;
ylim([ymin,ymax]);
title('$u$ flat: $\log_{10}$(relative error)','interpreter','latex','fontsize',16);

subplot(2,2,2);
contourf(Xcurve,Ycurve,log10(abs(Vexact-Vflat_inside)+eps),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
%caxis([-16,-12]);
axis equal;
ylim([ymin,ymax]);
title('$v$ flat: $\log_{10}$(absolute error)','interpreter','latex','fontsize',16);

subplot(2,2,3);
contourf(Xcurve,Ycurve,log10(abs((Uexact-Ucurve)./max(max(abs(Uexact))))+eps),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
%caxis([-16,-12]);
axis equal;
ylim([ymin,ymax]);
title('$u$ curve: $\log_{10}$(relative error)','interpreter','latex','fontsize',16);

subplot(2,2,4);
contourf(Xcurve,Ycurve,log10(abs(Vexact-Vcurve)+eps),20);
colorbar('TickLabelInterpreter','latex','fontsize',16);
%caxis([-16,-12]);
axis equal;
ylim([ymin,ymax]);
title('$v$ curve: $\log_{10}$(absolute error)','interpreter','latex','fontsize',16);
%saveas(gca,'/home/david/Documents/Research/kernels/slip_tests/pipe/vel_err_amp_1e-1_alpha_1e-1_eta_inf_sq_matlab','epsc');

%% extra functions
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

function [solution,K] = solve_direct(problem,rhs)
N = length(problem.domain.z);
K = zeros(2*N+2,2*N+2);
XX = eye(2*N+2);

parfor i = 1:2*N+2
    K(:,i) = matvec_combined_slip(XX(:,i),problem);
end

X = K\rhs;
solution.q = [X(1:N), X(N+1:2*N)];
solution.u_avg = [X(end-1); X(end)];
solution.trim = 1;
solution.local_indices = 1:length(solution.q);
solution.problem = problem;
end
