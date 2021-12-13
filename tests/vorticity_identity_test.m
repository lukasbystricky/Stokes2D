%% Vorticity identity test
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

% set target points
x = linspace(-2.5, 2.5, 50);
y = linspace(1, -1, 50);

%% Single-layer vorticity
problem.eta = inf;  % just SLP 
solution.problem = problem;
solution.trim = 0;
n = -1i*problem.domain.zp./abs(problem.domain.zp);
q = 1i*n;
solution.q = [real(q), imag(q)];
omegaslp = evaluate_vorticity(solution, x, y, 'verbose', 0);

% on-surface evalution
%solution.local_indices = 1:length(solution.q);
%omegaslp_on = evaluate_vorticity_on_surface(solution, solution, 'surface');

%% Double-layer vorticity
problem.eta = 0;    % just DLP
solution.problem = problem;
q = 1i*problem.domain.z;
solution.q = [real(q), imag(q)];
solution.trim = 0;

omegadlp = evaluate_vorticity(solution, x, y, 'verbose', 0);

% % on-surface evalution
% solution.local_indices = 1:length(solution.q);
% omegadlp_on = evaluate_vorticity_on_surface(solution, solution, 'surface');

%% Calculate error
e_slp = zeros(size(omegaslp));
e_dlp = zeros(size(omegadlp));
for i = 1:length(omegaslp)
    e_slp(i) = min(abs(omegaslp(i)), abs(omegaslp(i) - 1));
    e_dlp(i) = min(abs(omegadlp(i)), abs(omegadlp(i) - 2));
end

% % on-surface
% e_slp_on = zeros(size(omegaslp_on));
% e_dlp_on = zeros(size(omegadlp_on));
% for i = 1:length(omegaslp_on)
%     e_slp_on(i) = min(abs(omegaslp_on(i)), abs(omegaslp_on(i) - 0.5));
%     e_dlp_on(i) = min(abs(omegadlp_on(i)), abs(omegadlp_on(i) - 1));
% end

%% Plot
% Single-layer vorticity
figure;
semilogy(x,abs(e_slp));
ylabel('error in single-layer vorticity identity');
set(gca,'xticklabel',{[]});

% on-surface
% figure;
% semilogy(solution.problem.domain.theta,abs(e_slp_on));
% xlabel('t');
% ylabel('surface limiting value single-layer vorticity identity');
% xlim([0,2*pi]);

% Double-layer vorticity
figure;
semilogy(x,abs(e_dlp));
ylabel('error in double-layer vorticity identity');
set(gca,'xticklabel',{[]});

% on-surface
% figure;
% semilogy(solution.problem.domain.theta,abs(e_dlp_on));
% xlabel('t');
% ylabel('surface limiting value double-layer vorticity identity');
% xlim([0,2*pi]);
