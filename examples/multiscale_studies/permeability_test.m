% Investigation of the decay of the velocity boundary layer for flow over a
% porous bed.

close all
clearvars
clc

% create input structure
input_params = default_input_params('darcy_study', 1);

% modify structure as needed, or add additional problem-dependent params
Lx = 1;
Ly = 1;
circles = false;

input_params.box_size = [Lx,Ly];
input_params.panels = 40;
input_params.plot_domain = 1;

% set up radii and centers, centers given as complex numbers (x+iy)
c = 0.4;% concentration
x_centers = 0;
y_centers = 0;
centers = x_centers(:) + 1i*y_centers(:);
input_params.centers = centers;

if circles
    input_params.radii = Lx*sqrt((1-c)/pi);    

    problem = circles_periodic(input_params);
else
    %prescribe semi-major axis a
    input_params.a = 0.4*Lx;
    input_params.b = Lx^2*c./pi;
    input_params.angles = 5*pi/12;
    
    problem = ellipses_periodic(input_params);
end

%% compute permeability for a single reference obstacle

K = zeros(2,2);

% solve for K_{11} and K_{12} by imposing pressure gradient in x direction
problem.pressure_gradient_x = 1/Lx;
problem.pressure_gradient_y = 0;

% solve the problem
solution = solve_stokes(problem);

% average velocity is computed already!
K(:,1) = solution.u_avg;

% solve for K_{21} and K_{22} by imposing pressure gradient in y direction
problem.pressure_gradient_x = 0;
problem.pressure_gradient_y = 1/Ly;

% solve the problem
solution = solve_stokes(problem);

K(:,2) = solution.u_avg;

%% compute permeability for layer in y direction

input_params.centers = [-1i*Lx/2; 1i*Lx/2];
input_params.box_size = [Lx,2*Ly];

if circles
    input_params.radii = Lx*sqrt((1-c)/pi)*ones(2,1);    

    problem = circles_periodic(input_params);
else
    %prescribe semi-major axis a
    input_params.a = 0.4*Lx*ones(2,1);
    input_params.b = Lx^2*c./pi*ones(2,1);
    input_params.angles = 5*pi/12*ones(2,1);
    
    problem = ellipses_periodic(input_params);
end

Kvert = zeros(2,2);

% solve for K_{11} and K_{12} by imposing pressure gradient in x direction
problem.pressure_gradient_x = 1/Lx;
problem.pressure_gradient_y = 0;

% solve the problem
solution = solve_stokes(problem);

% average velocity is computed already!
Kvert(:,1) = solution.u_avg;

% solve for K_{21} and K_{22} by imposing pressure gradient in y direction
problem.pressure_gradient_x = 0;
problem.pressure_gradient_y = 1/Lx;

% solve the problem
solution = solve_stokes(problem);

Kvert(:,2) = solution.u_avg;

%% compute average velocity/pressure gradient in each layer and compare to
%expected Darcy velocity

input_params.box_size = [Lx,Ly];

x_centers = 0;
y_centers = 0;
centers = x_centers(:) + 1i*y_centers(:);
input_params.centers = centers;

if circles
    input_params.radii = Lx*sqrt((1-c)/pi)*ones(1,1);    

    problem = circles_periodic(input_params);
else
    %prescribe semi-major axis a
    input_params.a = 0.4*Lx*ones(1,1);
    input_params.b = Lx^2*c./pi*ones(1,1);
    input_params.angles = 5*pi/12*ones(1,1);
    
    problem = ellipses_periodic(input_params);
end

problem.pressure_gradient_x = 1;
problem.pressure_gradient_y = 0;

solution = solve_stokes(problem);

[u_avg, p_avg, p_grad_avg] = compute_cell_averages(solution, 1, 1, 0);
u_avg = u_avg(:,1) + 1i*u_avg(:,2);

% compute expected Darcy velocity
u_expected = zeros(size(u_avg,1),2);
for i = 1:size(u_expected,1)
    u_expected(i,:) = (Kvert*p_grad_avg(i,:)')';
end

u_expected = u_expected(:,1) + 1i*u_expected(:,2);