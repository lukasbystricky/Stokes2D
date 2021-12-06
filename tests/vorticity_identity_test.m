% %% Vorticity double-layer identity test
% close all
% clearvars
% clc
% 
% % create input structure
% input_params = default_input_params('non-periodic starfish', 1);
% 
% % modify structure as needed
% input_params.box_size = [5,3];
% input_params.panels = 40;
% input_params.centers = 0;
% input_params.radii = 1;
% 
% problem = starfish_periodic(input_params);
% problem.pressure_gradient_x = 0;
% problem.pressure_gradient_y = 0;
% problem.periodic = 1;
% solution.forces = 0;    % just DLP
% solution.trim = 0;
% 
% % set target points
% x = linspace(-2.5, 2.5, 50);
% y = linspace(1, -1, 50);

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

%% Double-layer vorticity
problem.eta = 0;    % just DLP
solution.problem = problem;
q = 1i*problem.domain.z;
%q = problem.domain.z;
solution.q = [real(q), imag(q)];
solution.trim = 0;
%omegadlp = evaluate_vorticity(solution, x, y, 'verbose', 0)
[omegadlp,omegadlp_no_sq,~,~] = evaluate_vorticity(solution, x, y, 'verbose', 0)

% on-surface evalution
solution.local_indices = 1:length(solution.q);
%omegadlp_on = evaluate_vorticity_on_surface(solution, solution, 'surface');
omegadlp_on = 0;

%% Calculate error
e_dlp = zeros(size(omegadlp));
for i = 1:length(omegadlp)
    e_dlp(i) = min(abs(omegadlp(i)), abs(omegadlp(i) - 2));
end

% on-surface
e_dlp_on = zeros(size(omegadlp_on));
for i = 1:length(omegadlp_on)
    e_dlp_on(i) = min(abs(omegadlp_on(i)), abs(omegadlp_on(i) - 1));
end

%% Plot
% Double-layer pressure
figure;
semilogy(x,abs(e_dlp));
ylabel('error in double-layer vorticity identity');
set(gca,'xticklabel',{[]});

figure;
semilogy(solution.problem.domain.theta,abs(e_dlp_on));
xlabel('t');
ylabel('surface limiting value double-layer vorticity identity');
xlim([0,2*pi]);

figure;
hold on;
plot(problem.domain.z,'b');
plot(x,y,'r.');
axis equal;
grid on;

figure;
subplot(1,2,1);
semilogy(x,abs(e_dlp),'bx--');
ylabel('error in double-layer vorticity identity');
%set(gca,'xticklabel',{[]});
xlabel('target x coordinate');
grid on;
%axis square;

subplot(1,2,2);
plot(x,omegadlp,'bx--');
ylabel('double-layer vorticity (with special quadrature)');
%set(gca,'xticklabel',{[]});
xlabel('target x coordinate');
grid on;
%axis square;