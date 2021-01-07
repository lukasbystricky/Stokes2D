% Investigation of the decay of the velocity boundary layer for flow over a
% porous bed.

close all
clearvars
clc

% create input structure
input_params = default_input_params('darcy_study', 1);

% modify structure as needed, or add additional problem-dependent params
n_layers = 20;
Lx = 1;
Ly = Lx*n_layers;
circles = false;

input_params.box_size = [Lx,Ly];
input_params.panels = 20;
input_params.plot_domain = 1;
input_params.pressure_drop_x = 1;
input_params.pressure_drop_y = 0;

% set up radii and centers, centers given as complex numbers (x+iy)
c = 0.2;% concentration
x_centers = zeros(n_layers/2,1);
y_centers = (-(n_layers/4)*Lx:Lx:(n_layers/4-1)*Lx) + Lx/2;
centers = x_centers(:) + 1i*y_centers(:);
input_params.centers = centers(1:end);

if circles
    radii = Lx*sqrt(c/pi);    
    input_params.radii = radii*ones(length(input_params.centers),1);
    
    problem_full = circles_periodic(input_params);
else
    %prescribe semi-major axis a
    input_params.a = 0.4*Lx*ones(length(input_params.centers),1);
    input_params.b = Lx^2*c./(pi * input_params.a);
    input_params.angles = 1*pi/12*ones(length(input_params.centers),1);
    
    problem_full = ellipses_periodic(input_params);
end

% solve the full problem
solution_full = solve_stokes(problem_full);

%% compute permeability for a single reference obstacle
input_params.box_size = [Lx,Lx];
input_params.centers = centers(1) + 1i*Lx/2;
input_params.plot_domain = 0;

if circles
    input_params.radii = input_params.radii(1);
    
    problem = circles_periodic(input_params);
else
    %prescribe semi-major axis a
    input_params.a = input_params.a(1);
    input_params.b = input_params.b(1);
    input_params.angles = input_params.angles(1);
    
    problem = ellipses_periodic(input_params);
end

K = zeros(2,2);

% solve for K_{11} and K_{12} by imposing pressure gradient in x direction
problem.pressure_gradient_x = 1/Lx;
problem.pressure_gradient_y = 0;

% solve the problem
solution_tmp = solve_stokes(problem);

% average velocity is computed already!
K(:,1) = solution_tmp.u_avg;

% solve for K_{21} and K_{22} by imposing pressure gradient in y direction
problem.pressure_gradient_x = 0;
problem.pressure_gradient_y = 1/Lx;

% solve the problem
solution_tmp = solve_stokes(problem);

K(:,2) = solution_tmp.u_avg;

%% compute averages in each layer
bodies = 1;
[u_avg, p_avg, p_grad_avg, u_grad_avg] = compute_cell_averages(solution_full, 0, 0, Lx, Lx, bodies);
u_avg = u_avg(:,1,:) + 1i*u_avg(:,2,:);

alpha = zeros(size(offsets));
interface_layer = n_layers/4;
for i = 1:length(offsets)
    alpha(i) = real(u_avg(interface_layer+1,i) - u_avg(interface_layer,i))/u_grad_avg(interface_layer,2,1,i);
end
% p_avg(1:end/2) = p_avg(1:end/2);
% p_grad_avg(1:end/2,:) = p_grad_avg(1:end/2,:);

%% plot Stokes solution
Nx = 50;
Ny = Nx*n_layers;
x = linspace(-Lx/2, Lx/2, Nx);
y = linspace(n_layers/2, -n_layers/2, Ny);
[X, Y] = meshgrid(x, y);

[U, V, X, Y] = evaluate_velocity(solution_full, X, Y);
P = evaluate_pressure(solution_full, X, Y);
[Px, Py] = evaluate_pressure_gradient(solution_full, X, Y);

Xstokes = X;
Ystokes = Y;

figure()
subplot(1,5,1);
contourf(Xstokes, Ystokes, U);
axis equal
title('U');
colorbar;

subplot(1,5,2);
contourf(Xstokes, Ystokes, V);
axis equal
title('V');
colorbar;

subplot(1,5,3);
contourf(Xstokes, Ystokes, P);
axis equal
title('P');
colorbar;

subplot(1,5,4);
contourf(Xstokes, Ystokes, Px);
axis equal
title('P_x');
colorbar;

subplot(1,5,5);
contourf(Xstokes, Ystokes, Py);
axis equal
title('P_y');
colorbar;

%% plot homogenized solution

Uavg = zeros(size(U));
Vavg = zeros(size(V));
Pavg = zeros(size(P));
Pxavg = zeros(size(Px));
Pyavg = zeros(size(Py));

for i = 1:length(u_avg)
    indices_tmp = (i-1)*Nx+1:i*Nx;
    
    Uavg(indices_tmp,:) = real(u_avg(end-i+1));
    Vavg(indices_tmp,:) = imag(u_avg(end-i+1));
    Pavg(indices_tmp,:) = p_avg(end-i+1);
    Pxavg(indices_tmp,:) = p_grad_avg(end-i+1,1);
    Pyavg(indices_tmp,:) = p_grad_avg(end-i+1,2);
end

Xhomogenized = X;
Yhomogenized = Y;

figure()

subplot(1,5,1);
contourf(Xhomogenized, Yhomogenized, Uavg);
axis equal
title('U^d');
colorbar;

subplot(1,5,2);
contourf(Xhomogenized, Yhomogenized, Vavg);
axis equal
title('V^d');
colorbar;

subplot(1,5,3);
contourf(Xhomogenized, Yhomogenized, Pavg);
axis equal
title('P^d');
colorbar;

subplot(1,5,4);
contourf(Xhomogenized, Yhomogenized, Pxavg);
axis equal
title('P_x^d');
colorbar;

subplot(1,5,5);
contourf(Xhomogenized, Yhomogenized, Pyavg);
axis equal
title('P_y^d');
colorbar;
