% Investigation of the decay of the velocity boundary layer for flow over a
% porous bed.

close all
clearvars
clc

% create input structure
input_params = default_input_params('darcy_study', 1);

% modify structure as needed, or add additional problem-dependent params
n_layers = 10;
Lx = 1;
Ly = Lx*n_layers;
circles = true;
array = false;

input_params.box_size = [Lx,Ly];
input_params.panels = 20;
input_params.plot_domain = 1;
input_params.pressure_drop_x = 0;
input_params.pressure_drop_y = n_layers;

% set up radii and centers, centers given as complex numbers (x+iy)
c = 0.4;% concentration
x_centers = zeros(n_layers,1);
y_centers = -n_layers/2*Lx:Lx:(n_layers-1)*Lx/2;
centers = x_centers(:) + 1i*(y_centers(:)+Lx/2);
input_params.centers = centers(1:end/2);

if circles
    radii = Lx*sqrt(c/pi);
    ell = Lx;
    
    if array
        radii = radii/2;
        ell = 0.5*Lx;
        
        x_centers = -0.25*ones(2*n_layers,1);
        y_centers = [y_centers + 0.25; y_centers - 0.25];
        
        centers = x_centers(:) + 1i*(y_centers(:)+Lx/2);
        
        input_params.centers = [centers(1:end/2); centers(1:end/2)+0.5];
    end

    
    input_params.radii = radii*ones(length(input_params.centers),1);
    
    problem_full = circles_periodic(input_params);
else
    %prescribe semi-major axis a
    input_params.a = 0.4*Lx*ones(length(input_params.centers),1);
    input_params.b = Lx^2*c./(pi * input_params.a);
    input_params.angles = 4*pi/12*ones(length(input_params.centers),1);
    
    problem_full = ellipses_periodic(input_params);
end

% solve the full problem
solution_full = solve_stokes(problem_full);

%% compute permeability for a single reference obstacle
input_params.box_size = [Lx,Lx];
input_params.centers = 0;
input_params.plot_domain = 1;

if circles
    
    if array
        input_params.radii = 2*input_params.radii(1);
    else        
        input_params.radii = input_params.radii(1);
    end
    
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

% scale by length scales
K = K*(ell/Lx)^2;

%% compute averages in each layer
[u_avg, p_avg, p_grad_avg, u_laplace_avg, x_cen, y_cen] = compute_cell_averages(solution_full, 1, Ly, 0);
u_avg = u_avg(:,1,:) + 1i*u_avg(:,2,:);

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

%% plot averaged solution
n_averages = size(y_cen,2);

figure();

for i = 1:n_averages
    subplot(5,1,1)
    hold on
    stairs(y_cen(:,i), real(u_avg(:,:,i)));
    title('tangential velocity');
    
    subplot(5,1,2)
    hold on
    stairs(y_cen(:,i), imag(u_avg(:,:,i)));
    title('normal velocity');
    
    subplot(5,1,3)
    hold on
    stairs(y_cen(:,i), p_avg(:,i));
    title('pressure');
    
    subplot(5,1,4)
    hold on
    stairs(y_cen(:,i), p_grad_avg(:,1,i));
    title('tangential pressure gradient');
    
    subplot(5,1,5)
    hold on
    stairs(y_cen(:,i), p_grad_avg(:,2,i));
    title('normal pressure gradient');
end
%% plot homogenized solution
% 
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

%% compute expected Darcy velocity 
u_expected = zeros(size(u_avg,1),2);
p_expected = zeros(size(u_avg,1),1);
px_expected = zeros(size(u_avg,1),1);
py_expected = zeros(size(u_avg,1),1);

for i = 1:n_layers/2
    u_expected(i,:) = (K*p_grad_avg(i,:)')';
    p_expected(i) = nan;
    px_expected(i) = nan;
    py_expected(i) = nan;
end

% expected velocity above interface is the velocity evaluated at the cell
% center
for i = 1:n_layers/2
    [utmp, vtmp] = evaluate_velocity(solution_full, 0, (i-1)+0.5);
    [pxtmp,pytmp] = evaluate_pressure_gradient(solution_full, 0, (i-1)+0.5);
    ptmp = evaluate_pressure(solution_full, 0, (i-1)+0.5);
    
    u_expected(n_layers/2 + i,1) = utmp;
    u_expected(n_layers/2 + i,2) = vtmp;
    p_expected(n_layers/2 + i) = ptmp;
    px_expected(n_layers/2 +i) = pxtmp;
    py_expected(n_layers/2 +i) = pytmp;
end

u_expected = u_expected(:,1) + 1i*u_expected(:,2);

figure();
subplot(5,1,1);
plot(real(u_avg));
hold on
plot(real(u_expected));
title('velocity (horizontal)');
legend('average velocity', 'expected velocity');

subplot(5,1,2);
plot(imag(u_avg));
hold on
plot(imag(u_expected));
title('velocity (vertical)');
legend('average velocity', 'expected velocity');

subplot(5,1,3);
plot(log10(abs(u_avg - u_expected)));
title('error (log_{10})');

subplot(5,1,4);
plot(p_avg);
hold on
plot(p_expected);
title('pressure');

subplot(5,1,5);
plot(p_grad_avg(:,2));
hold on
plot(px_expected);
title('pressure gradient (vertical)');
% 
% 
% %% compute force on interface
% h = Lx/200;
% x = h:h:Lx;
% y = 0*ones(size(x));
% z = x + 1i*y;
% zp = ones(size(x))/Lx;
% w = h*ones(size(x));
% 
% f = compute_force(solution_full, z, zp, w);
