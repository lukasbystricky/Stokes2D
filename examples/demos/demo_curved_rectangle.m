close all
clearvars
clc

% create input structure
input_params = default_input_params('curved_rectangle', 0);

% modify structure as needed, or add additional problem-dependent params
input_params.plot_domain = 0;

% position of the four corners
input_params.corners = [-0.5+1i*0.5; 0.5+1i*0.5; 0.5-1i*0.5; -0.5-1i*0.5];

% boundary condition for the walls ordered as (top, right, bottom, left).
% Here it is assumed that the walls are disretzied with the same number of
% panels
input_params.boundary_conditions = @(z) [0*ones(length(z)/4,1);
                                         0*ones(length(z)/4,1);
                                         -1*ones(length(z)/4,1);
                                         0*ones(length(z)/4,1)];
                                     
% number of panels on the walls (can be specified individually for each
% wall as e.g. [4 5 3 10])
input_params.panels = 5;

% naively subdivide panels closest to the corners (set 0 for no divisions)
input_params.nsub = 5;

% create problem
problem = curved_rectangle(input_params);

%% solve problem
solution = solve_stokes(problem,'fmm',0);

% trim evaluation if point is outside the domain
solution.trim = 1;

%% evaluate velocity inside domain
[U,V,X,Y] = evaluate_velocity(solution,200,'fmm',1,'verbose',0);

%% evaluate velocity along a line inside the domain
npan = 2;
line{1} = @(T) geometry_periodic_channel(@(t) 0*ones(size(t)), ...
    @(t) zeros(size(t)), @(t) zeros(size(t)),T, 1, 1);
line_domain = discretize_domain(line,npan);
x = real(line_domain.z);
y = imag(line_domain.z);

[Uline,Vline] = evaluate_velocity(solution,x,y,'fmm',1,'verbose',0);

%% integrate quantity along the same line
Uint = sum(Uline.*line_domain.wazp)
Vint = sum(Vline.*line_domain.wazp)

%% plot
plot_domain(problem);
hold on;
plot(x,y,'mx');
xlabel('x');
ylabel('y');
grid on;
legend('collocation points','normal vectors','line integral points');

figure;
subplot(1,2,1);
mesh(X,Y,U,'facecolor','interp');
view(2);
colorbar;
xlabel('x');
ylabel('y');
title('U');
grid off;
axis equal;

subplot(1,2,2);
mesh(X,Y,V,'facecolor','interp');
view(2);
colorbar;
xlabel('x');
ylabel('y');
title('V');
grid off;
axis equal;

figure;
quiver(X,Y,U,V,3);
xlabel('x');
ylabel('y');
title('velocity components');
axis equal;

figure;
hold on;
plot(x,Uline);
plot(x,Vline);
legend('U','V');
xlabel('x');
title('velocity along horizontal line at y=0');
grid on;
