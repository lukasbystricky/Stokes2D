%% Run stresslet ID tests
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

problem = starfish_periodic(input_params);
problem.pressure_gradient_x = 0;
problem.pressure_gradient_y = 0;

% set target points
x = linspace(-2.5, 2.5, 50);
y = linspace(1, -1, 50);

%% Single-layer pressure
% Here we set the density function to be the normal vector
problem.eta = inf; % just SLP 
solution.problem = problem;
solution.trim = 0;
q = 1i*problem.domain.zp./abs(problem.domain.zp);
solution.q = [real(q), imag(q)];
Pslp = evaluate_pressure(solution, x, y, 'verbose', 0);

% on-surface evalution
solution.local_indices = 1:length(solution.q);
Pslp_on = evaluate_pressure_on_surface(solution, solution, 'surface');

%% Double-layer pressure
solution.problem.eta = 0; % just DLP
q = problem.domain.z;
solution.q = [real(q), imag(q)];
Pdlp = evaluate_pressure(solution, x, y, 'verbose', 0);

% on-surface evaluation
solution.local_indices = 1:length(solution.q);
Pdlp_on = evaluate_pressure_on_surface(solution, solution, 'surface');

%% Calculate error
e_slp = zeros(size(Pslp));
e_dlp = zeros(size(Pdlp));
for i = 1:length(Pslp)
    e_slp(i) = min(Pslp(i), abs(Pslp(i) - 1));
    e_dlp(i) = min(abs(Pdlp(i)), abs(Pdlp(i) - 2));
end

% on-surface
e_slp_on = zeros(size(Pslp_on));
e_dlp_on = zeros(size(Pdlp_on));
for i = 1:length(Pslp_on)
    e_slp_on(i) = min(Pslp_on(i), abs(Pslp_on(i) - 0.5));
    e_dlp_on(i) = min(abs(Pdlp_on(i)), abs(Pdlp_on(i) - 1));
end

%% Plot
% Single-layer pressure
figure;
semilogy(x,abs(e_slp));
ylabel('error in single-layer pressure identity');
set(gca,'xticklabel',{[]});

figure;
semilogy(solution.problem.domain.theta,abs(e_slp_on));
xlabel('t');
ylabel('surface limiting value single-layer pressure identity');
xlim([0,2*pi]);

% Double-layer pressure
figure;
semilogy(x,abs(e_dlp));
ylabel('error in double-layer pressure identity');
set(gca,'xticklabel',{[]});

figure;
semilogy(solution.problem.domain.theta,abs(e_dlp_on));
xlabel('t');
ylabel('surface limiting value double-layer pressure identity');
xlim([0,2*pi]);
