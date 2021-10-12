close all
clearvars
clc

% create input structure
input_params = default_input_params('pouseuille_demo', 1);

% modify structure as needed
input_params.box_size = [5,5];
input_params.h = 0.5;
input_params.panels = 40;
input_params.eta = 1;
input_params.plot_domain = 0;

problem = flat_pipe_periodic(input_params);

%% solve the problem
solution = solve_stokes(problem,'fmm',0);
solution.local_indices = 1:length(solution.q);

%% evaulate on surface
Tsurf = problem.domain.theta;
[Usurf, Vsurf] = evaluate_velocity_on_surface(solution, solution);
[Uxsurf, Uysurf, Vxsurf, Vysurf] = evaluate_velocity_gradient_on_surface(solution, solution, 'fluid');
Psurf = evaluate_pressure_on_surface(solution, solution, 'fluid');
[Pxsurf, Pysurf] = evaluate_pressure_gradient_on_surface(solution, solution, 'fluid');

%% evaluate near surface
%Tnear = linspace(0,2*pi,600)';
Tnear = Tsurf;
ztar = problem.domain.z - 1e-3*1i;
%ztar = [problem.walls{1}(Tnear); problem.walls{2}(Tnear)];
%Tnear = linspace(0,2*pi,2*600)';
xtar = real(ztar);
ytar = imag(ztar);

[Unear, Vnear, ~, ~, ~, ~] = evaluate_velocity(solution, xtar, ytar, 'fmm', 0, 'verbose', 0);
[Uxnear, Uynear, Vxnear, Vynear, ~, ~, ~, ~] = evaluate_velocity_gradient(solution, xtar, ytar);
[Pnear, ~, ~, ~] = evaluate_pressure(solution, xtar, ytar, 'fmm', 0, 'verbose', 0);
[Pxnear, Pynear] = evaluate_pressure_gradient(solution, xtar, ytar);

%% plot
close all;
figure;
subplot(1,2,1);
plot(Tsurf,Usurf,'x');
hold on;
plot(Tnear,Unear,'o');
legend('On surface','Near surface');
title('U');

subplot(1,2,2);
plot(Tsurf,Vsurf,'x');
hold on;
plot(Tnear,Vnear,'o');
legend('On surface','Near surface');
title('V');

figure;
subplot(2,2,1);
plot(Tsurf,Uxsurf,'x');
hold on;
plot(Tnear,Uxnear,'o');
legend('On surface','Near surface');
title('Ux');

subplot(2,2,2);
plot(Tsurf,Uysurf,'x');
hold on;
plot(Tnear,Uynear,'o');
legend('On surface','Near surface');
title('Uy');

subplot(2,2,3);
plot(Tsurf,Vxsurf,'x');
hold on;
plot(Tnear,Vxnear,'o');
legend('On surface','Near surface');
title('Vx');

subplot(2,2,4);
plot(Tsurf,Vysurf,'x');
hold on;
plot(Tnear,Vynear,'o');
legend('On surface','Near surface');
title('Vy');

figure;
plot(Tsurf,Psurf,'x');
hold on;
plot(Tnear,Pnear,'o');
legend('On surface','Near surface');
title('P');

figure;
subplot(1,2,1);
plot(Tsurf,Pxsurf,'x');
hold on;
plot(Tnear,Pxnear,'o');
legend('On surface','Near surface');
title('Px');

subplot(1,2,2);
plot(Tsurf,Pysurf,'x');
hold on;
plot(Tnear,Pynear,'o');
legend('On surface','Near surface');
title('Py');