%% Run stresslet ID tests
%close all
clearvars
clc

% create input structure
input_params = default_input_params('periodic_circles', 1);

% modify structure as needed
input_params.box_size = [5,3];
input_params.panels = 40;
input_params.centers = 0;
input_params.radii = 1;

problem =  starfish_periodic(input_params);
%problem = circles_periodic(input_params);
problem.pressure_gradient_x = 0;
problem.pressure_gradient_y = 0;

% set target points
x = linspace(-2.5, 2.5, 50);
y = linspace(1, -1, 50);

%% Single-layer pressure
% Here we set the density function to be the normal vector

problem.eta = inf; % just SLP 
solution.problem = problem;
q = 1i*problem.domain.zp./abs(problem.domain.zp);
solution.q = [real(q), imag(q)];
Pslp = evaluate_pressure(solution, x, y);

% on-surface evalution
solution.local_indices = 1:length(solution.q);
Pslp_on = evaluate_pressure_on_surface(solution, solution, 'surface');

% figure(1);
% semilogy(solution.problem.domain.theta, abs(Pslp_on));
% hold on
solution.problem.eta = 0; % just DLP
q = problem.domain.z;
solution.q = [real(q), imag(q)];
Pdlp = evaluate_pressure(solution, x, y);

% on-surface evalution
solution.local_indices = 1:length(solution.q);
Pdlp_on = evaluate_pressure_on_surface(solution, solution, 'surface');

e_slp = zeros(size(Pslp));
e_dlp = zeros(size(Pslp));
% calculate error
for i = 1:length(Pslp)
    e_slp(i) = min(Pslp(i), abs(Pslp(i) - 1));
    e_dlp(i) = min(abs(Pdlp(i)), abs(Pdlp(i) - 1));
end

addpath('../matlab2tikz/src');
% plot(solution.problem.domain.z);
% hold on
% plot(x,y);
% axis off
% axis equal
% matlab2tikz('standalone', true, 'domain_evaluation.tex');

figure(1);
semilogy(x,abs(e_slp));
hold on
set(gca,'xticklabel',{[]});
%matlab2tikz('standalone', true, 'pressure_slp_identity.tex');

figure(2);
semilogy(x,abs(e_dlp));
hold on
set(gca,'xticklabel',{[]});

%figure(2)
%semilogy(solution.problem.domain.theta, abs(Pdlp_on));
%hold on
%figure();
%plot(x, e_dlp);
% 
% figure()
% plot(Pdlp_on);


