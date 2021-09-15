% Investigation of the decay of the velocity boundary layer for flow over a
% porous bed.

function [alpha, alpha_saffman, alpha_bar, alpha_saffman_bar, U, K, interface,...
                infiltration_angle] = compute_alpha(pressure_angle,...
            interface_offset, volume_fraction, theta)

% close all
% clearvars
% clc

% create input structure
input_params = default_input_params('darcy_study', 1);

% modify structure as needed, or add additional problem-dependent params
n_layers_free = 20;
n_pores = 20;
Lx = 1;
Ly = Lx*(n_layers_free + n_pores);
circles = 0;

input_params.box_size = [Lx,Ly];
input_params.panels = 20;
input_params.plot_domain = 1;
input_params.eta = 1;

% set pressure drop based on angle
input_params.pressure_drop_x = sqrt(1/(1+tan(pressure_angle)^2))*Lx;
input_params.pressure_drop_y = tan(pressure_angle)*input_params.pressure_drop_x*Ly;


% set up radii and centers, centers given as complex numbers (x+iy)
c = volume_fraction;% concentration
x_centers = zeros(n_pores,1);
y_centers = (-(n_pores/2)*Lx:Lx:(n_pores/2-1)*Lx) + Lx/2;
centers = x_centers(:) + 1i*y_centers(:);
input_params.centers = centers(1:end);

if circles
    radii = Lx*sqrt(c/pi);    
    input_params.radii = radii*ones(length(input_params.centers),1);
    
    problem_full = circles_periodic(input_params);
else
    %prescribe semi-major axis a
    input_params.b = Lx * sqrt(c /(2*pi))*ones(length(input_params.centers),1);
    input_params.a = 2*input_params.b;
    input_params.angles = theta*ones(length(input_params.centers),1);
    
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
problem.pressure_gradient_x = 1;
problem.pressure_gradient_y = 0;

% solve the problem
solution_tmp = solve_stokes(problem);

% average velocity is computed already!
K(:,1) = solution_tmp.u_avg;

% solve for K_{21} and K_{22} by imposing pressure gradient in y direction
problem.pressure_gradient_x = 0;
problem.pressure_gradient_y = 1;

% solve the problem
solution_tmp = solve_stokes(problem);

K(:,2) = solution_tmp.u_avg;

%% compute averages in each layer

n_boundary_layers = 0;
interface = max(imag(problem_full.domain.z)) + interface_offset;
%interface = 10;
y_darcy = interface - Lx/2 - n_boundary_layers*Lx : Lx :  Ly/2;
[u_avg, ~, p_grad_avg, u_grad_avg] = compute_cell_averages(solution_full, ...
                0, y_darcy, Lx, Lx);
            
u_avg = u_avg(:,1) + 1i*u_avg(:,2);
u_expected = zeros(size(u_avg));

u_stokes = u_avg(2);
u_darcy = u_avg(1);
uy_stokes = u_grad_avg(2,2,1);

alpha = uy_stokes/real(u_stokes - u_darcy) * sqrt(trace(abs(K)));
alpha_saffman = uy_stokes/real(u_stokes) * sqrt(trace(abs(K)));

x = linspace(-Ly/2, Ly/2, 30);
[X,Y] = meshgrid(x, interface);

U = evaluate_velocity(solution_full, X, Y);
[~, Uy] = evaluate_velocity_gradient(solution_full, X, Y);
u_stokes1 = mean(U);
uy_stokes1 = mean(Uy);

alpha_bar = uy_stokes1/real(u_stokes1 - u_darcy) * sqrt(trace(abs(K)));
alpha_saffman_bar = uy_stokes1/real(u_stokes1) * sqrt(trace(abs(K)));

infiltration_angle = atan(imag(u_avg(end))/real(u_avg(end)));

Ny = 50;
x = 0;
y = linspace(interface, Ly/2, Ny);
[X, Y] = meshgrid(x, y);

[U,V] = evaluate_velocity(solution_full, X, Y);
[~, Uy] = evaluate_velocity_gradient(solution_full, X, Y);


% for i = 1:length(y_darcy)
%    utmp = K*p_grad_avg(i,:)';
%    u_expected(i) = utmp(1) + 1i*utmp(2);
% end

% figure();
% subplot(2,1,1)
% plot(log10(abs(real(u_avg))), 'b');
% hold on
% plot(log10(abs(real(u_expected))), '--b');
% 
% subplot(2,1,2);
% plot(log10(abs(real(u_expected - u_avg))));

%% plot Stokes solution
% Nx = 50;
% Ny = Nx*(n_layers_free + n_pores);
% x = linspace(-Lx/2, Lx/2, Nx);
% y = linspace(-Ly/2, Ly/2, Ny);
% [X, Y] = meshgrid(x, y);
% 
% [U, V, X, Y] = evaluate_velocity(solution_full, X, Y);
% P = evaluate_pressure(solution_full, X, Y);
% [Px, Py] = evaluate_pressure_gradient(solution_full, X, Y);
% 
% Xstokes = X;
% Ystokes = Y;
% 
% figure()
% subplot(1,5,1);
% contourf(Xstokes, Ystokes, U);
% axis equal
% title('U');
% colorbar;
% 
% subplot(1,5,2);
% contourf(Xstokes, Ystokes, V);
% axis equal
% title('V');
% colorbar;
% 
% subplot(1,5,3);
% contourf(Xstokes, Ystokes, P);
% axis equal
% title('P');
% colorbar;
% 
% subplot(1,5,4);
% contourf(Xstokes, Ystokes, Px);
% axis equal
% title('P_x');
% colorbar;
% 
% subplot(1,5,5);
% contourf(Xstokes, Ystokes, Py);
% axis equal
% title('P_y');
% colorbar;
% 
% %% plot homogenized solution
% 
% Uavg = zeros(size(U));
% Vavg = zeros(size(V));
% Pavg = zeros(size(P));
% Pxavg = zeros(size(Px));
% Pyavg = zeros(size(Py));
% 
% for i = 1:length(u_avg)
%     indices_tmp = (i-1)*Nx+1:i*Nx;
%     
%     Uavg(indices_tmp,:) = real(u_avg(end-i+1));
%     Vavg(indices_tmp,:) = imag(u_avg(end-i+1));
%     Pavg(indices_tmp,:) = p_avg(end-i+1);
%     Pxavg(indices_tmp,:) = p_grad_avg(end-i+1,1);
%     Pyavg(indices_tmp,:) = p_grad_avg(end-i+1,2);
% end
% 
% Xhomogenized = X;
% Yhomogenized = Y;
% 
% figure()
% 
% subplot(1,5,1);
% contourf(Xhomogenized, Yhomogenized, Uavg);
% axis equal
% title('U^d');
% colorbar;
% 
% subplot(1,5,2);
% contourf(Xhomogenized, Yhomogenized, Vavg);
% axis equal
% title('V^d');
% colorbar;
% 
% subplot(1,5,3);
% contourf(Xhomogenized, Yhomogenized, Pavg);
% axis equal
% title('P^d');
% colorbar;
% 
% subplot(1,5,4);
% contourf(Xhomogenized, Yhomogenized, Pxavg);
% axis equal
% title('P_x^d');
% colorbar;
% 
% subplot(1,5,5);
% contourf(Xhomogenized, Yhomogenized, Pyavg);
% axis equal
% title('P_y^d');
% colorbar;
