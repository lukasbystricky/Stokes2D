% DESCRIPTION

close all;
clearvars;
clc;
format long;
set(groot,'defaultAxesTickLabelInterpreter','latex');

% parameters
alpha = 0;
amplitude = 1e-1;

% create input structure
input_params = default_input_params('slip_pipe_test', 1);

% modify structure as needed
input_params.box_size = [2,2];
input_params.h = 0.5;
input_params.panels = 20;
input_params.eta = inf;
input_params.plot_domain = 0;
input_params.slip = 1;
input_params.d = 0.5;
input_params.alpha = alpha;
input_params.A = amplitude;

% initialize problem
problem = sine_pipe_periodic(input_params);
domain = problem.domain;

% normalized tangent
nu1 = real(domain.zp)./abs(domain.zp);
nu2 = imag(domain.zp)./abs(domain.zp);

% normalized normal
n1 = real(-1i*domain.zp)./abs(domain.zp);
n2 = imag(-1i*domain.zp)./abs(domain.zp);

%% solve problem
N = length(problem.domain.z);

rhs = [zeros(N/2,1);
       zeros(N/2,1);
       problem.pressure_gradient_x; problem.pressure_gradient_y];

solution = solve_stokes(problem,rhs,'fmm',0);
solution.local_indices = 1:length(solution.q);

%% on-surface
[u1,u2] = evaluate_velocity_on_surface(solution,solution);
[sxx,sxy,syx,syy] = evaluate_stress_on_surface(solution,solution,'fluid');

%% boundary condition
% normal direction
udotn = u1.*n1 + u2.*n2;

% tangential direction
udotnu = u1.*nu1 + u2.*nu2;
t1 = n1.*sxx + n2.*syx;
t2 = n1.*sxy + n2.*syy;
nudott = t1.*nu1 + t2.*nu2;

% final
g1 = udotn;
g2 = udotnu + alpha*nudott;

%% plot
figure;
subplot(1,2,1);
plot(g1);
title('g1');

subplot(1,2,2);
plot(g2);
title('g2');
