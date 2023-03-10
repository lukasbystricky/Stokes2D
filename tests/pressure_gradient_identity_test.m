%% Pressure gradient identity test
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

%% Single-layer pressure gradient
problem.eta = inf;  % just SLP 
solution.problem = problem;
solution.trim = 0;
n = -1i*problem.domain.zp./abs(problem.domain.zp);
q = n.*problem.domain.z;
solution.q = [real(q), imag(q)];
dpslp = evaluate_pressure_gradient(solution, x, y, 'verbose', 0);

% on-surface evalution
solution.local_indices = 1:length(solution.q);
dpslp_on = evaluate_pressure_gradient_on_surface(solution, solution, 'surface');

%% Double-layer pressure gradient
problem.eta = 0;    % just DLP
solution.problem = problem;
q = (problem.domain.z).^2;
solution.q = [real(q), imag(q)];
solution.trim = 0;
dpdlp = evaluate_pressure_gradient(solution, x, y, 'verbose', 0);

% on-surface evalution
solution.local_indices = 1:length(solution.q);
dpdlp_on = evaluate_pressure_gradient_on_surface(solution, solution, 'surface');

%% Calculate error
e_slp = zeros(size(dpslp));
e_dlp = zeros(size(dpdlp));
for i = 1:length(dpslp)
    e_slp(i) = min(abs(dpslp(i)), abs(dpslp(i) - -1));
    e_dlp(i) = min(abs(dpdlp(i)), abs(dpdlp(i) - 2));
end

% on-surface
e_slp_on = zeros(size(dpslp_on));
e_dlp_on = zeros(size(dpdlp_on));
for i = 1:length(dpslp_on)
    e_slp_on(i) = min(abs(dpslp_on(i)), abs(dpslp_on(i) - -0.5));
    e_dlp_on(i) = min(abs(dpdlp_on(i)), abs(dpdlp_on(i) - 1));
end

%% Plot
% Single-layer pressure gradient
figure;
semilogy(x,abs(e_slp));
ylabel('error in single-layer pressure gradient identity');
set(gca,'xticklabel',{[]});

% on-surface
figure;
semilogy(solution.problem.domain.theta,abs(e_slp_on));
xlabel('t');
ylabel('error in surface limiting value single-layer pressure gradient identity');
xlim([0,2*pi]);

% Double-layer pressure gradient
figure;
semilogy(x,abs(e_dlp));
ylabel('error in double-layer pressure gradient identity');
set(gca,'xticklabel',{[]});

% on-surface
figure;
semilogy(solution.problem.domain.theta,abs(e_dlp_on));
xlabel('t');
ylabel('error in surface limiting value double-layer pressure gradient identity');
xlim([0,2*pi]);
