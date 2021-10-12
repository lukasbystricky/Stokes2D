% Run stresslet ID tests
close all
clearvars
clc

% create input structure
input_params = default_input_params('periodic_starfish', 1);

% modify structure as needed
input_params.box_size = [5,3];
input_params.panels = 40;
input_params.centers = 0;
input_params.radii = 1;
input_params.plot_domain = 0;

problem = starfish_periodic(input_params);
problem.pressure_gradient_x = 0;
problem.pressure_gradient_y = 0;

%% Double-layer velocity
solution.problem = problem;
solution.trim = 0;      % do not remove values outside the domain
solution.problem.eta = 0;   % just DLP
solution.q = [ones(size(problem.domain.z)), zeros(size(problem.domain.z))];
solution.u_avg = [0; 0];

Lx = solution.problem.domain.Lx;
Ly = solution.problem.domain.Ly;
xsrc = real(solution.problem.domain.z);
ysrc = imag(solution.problem.domain.z);
n1 = real(-1i*solution.problem.domain.zp)./abs(solution.problem.domain.zp);
n2 = imag(-1i*solution.problem.domain.zp)./abs(solution.problem.domain.zp);
weights = solution.problem.domain.wazp;
domain = solution.problem.domain;

% set target points
M = 100;
x = linspace(-2.5, 2.5, M);
y = linspace(2.5, -2.5, M);
[X,Y] = meshgrid(x,y);
Xtar_sq = mod(X+Lx/2,Lx)-Lx/2;
Ytar_sq = mod(Y+Ly/2,Ly)-Ly/2;

[Uc, Vc, X, Y, U, V] = evaluate_velocity(solution, X, Y, 'fmm', 0, 'verbose', 0);

solution.local_indices = 1:length(solution.q);
Psurface = evaluate_pressure_on_surface(solution, solution, 'fluid');

%% Find points inside domain by applying stresslet identity
[test1, test2] = StokesDLP_ewald_2p(xsrc, ysrc, X(:), Y(:), n1, n2,...
    ones(length(n1),1).*weights, zeros(length(n1),1).*weights, Lx, Ly);

% correct using special quadrature
[test,~] = mex_SQ_dlp(Xtar_sq(:)+1i*(Ytar_sq(:)+1e-60), domain.z,...
                domain.zp, domain.quad_weights, ...
                domain.panel_breaks, domain.wazp, domain.z32, domain.zp32,...
                domain.quad_weights32, domain.wazp32,ones(length(n1),1) + 1e-14*1i,...
                test1 + 1i*test2,domain.mean_panel_length,domain.extra.gridSolidmat, ...
                domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
                domain.reference_cell,solution.problem.periodic);

% anything that is greater than 0 is inside the starfish
inside = find(solution.problem.stresslet_id_test(real(test)) == 1);
outside = setdiff(1:length(X(:)),inside)';

%% Separate velocity components outside and inside the domain
Uc_out = Uc(:);
Uc_out(inside) = nan;
Uc_in = Uc(:);
Uc_in(outside) = nan;

Vc_out = Vc(:);
Vc_out(inside) = nan;
Vc_in = Vc(:);
Vc_in(outside) = nan;

X_out = X;
X_out(inside) = nan;
X_in = X;
X_in(outside) = nan;

Y_out = Y;
Y_out(inside) = nan;
Y_in = Y;
Y_in(outside) = nan;

Uc_out = reshape(Uc_out,size(X));
Uc_in = reshape(Uc_in,size(X));

Vc_out = reshape(Vc_out,size(X));
Vc_in = reshape(Vc_in,size(X));

%% Plot
% u
figure;
subplot(1,3,1);
plot(problem.domain.z,'k');
hold on;
contourf(X_out,Y_out,log10(abs(Uc_out)+eps));
colorbar;
axis equal;
title('u (outside): log_{10} error');

subplot(1,3,2);
plot(problem.domain.z,'k');
hold on;
contourf(X_in,Y_in,log10(abs(Uc_in-1)+eps));
colorbar;
axis equal;
xlim([min(x),max(x)]);
ylim([min(y),max(y)]);
title('u (inside): log_{10} error');

subplot(1,3,3);
plot(problem.domain.z,'k');
hold on;
contourf(X_out,Y_out,log10(abs(Uc_out)+eps));
hold on;
contourf(X_in,Y_in,log10(abs(Uc_in-1)+eps));
colorbar;
axis equal;
title('u (outside + inside): log_{10} error');

% v
figure;
subplot(1,3,1);
plot(problem.domain.z,'k');
hold on;
contourf(X_out,Y_out,Vc_out);
colorbar;
axis equal;
title('v (outside)');

subplot(1,3,2);
plot(problem.domain.z,'k');
hold on;
contourf(X_in,Y_in,Vc_in);
colorbar;
axis equal;
xlim([min(x),max(x)]);
ylim([min(y),max(y)]);
title('v (inside)');

subplot(1,3,3);
plot(problem.domain.z,'k');
hold on;
contourf(X_out,Y_out,Vc_out);
hold on;
contourf(X_in,Y_in,Vc_in);
colorbar;
axis equal
title('v (outside + inside)');