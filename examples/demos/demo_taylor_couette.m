% Computes the velocity between two rotating cylinders, i.e. Taylor-Couette
% flow. The outer cylinder is fixed, while the inner cylinder rotates
% counter clockwise with angular velocity omega. Then, the exact solution
% is known, with angular velocity given by A*r + B/r, where
% A = -omega/(ro^2*(1/ri^2-1/ro^2)),
% B = omega/(1/ri^2-1/ro^2),
% and where ro and ri denote the radius of the outer and inner cylinder
% respectively
close all
clearvars
clc

test = 1;

% create input structure
input_params = default_input_params('couette', 0);

% modify structure as needed, or add additional problem-dependent params
input_params.panels = [10, 10];
input_params.plot_domain = 0;

% set up radii and centers, centers given as complex numbers (x+iy)
input_params.radii = [2, 0.5]; % the outer wall must be the first radii
input_params.centers = [0, 0];
input_params.omega = 1; % angular velocity of inner cylinder
input_params.slip = 0;

problem = couette(input_params);

% solve the problem
z = problem.domain.z;

rhs = [real(problem.boundary_conditions(z));... 
       imag(problem.boundary_conditions(z))];
solution = solve_stokes(problem,rhs,'fmm',1);
solution.local_indices = 1:length(solution.q);
solution.trim = 1;

%% off-surface
% grid in polar coordinates with M number of radial and angular points
M = 100;

[Uc, Vc, X, Y, U, V] = evaluate_velocity(solution, M, 'fmm', 1, 'verbose', 1);
[Pc, P, ~, ~] = evaluate_pressure(solution, X, Y, 'fmm', 0, 'verbose', 1);
[Uxc, Uyc, Vxc, Vyc, Ux, Uy, Vx, Vy] = evaluate_velocity_gradient(solution, X, Y);
[Pxc, Pyc, Px, Py, ~, ~] = evaluate_pressure_gradient(solution, X, Y);
[omegac, ~, ~, ~] = evaluate_vorticity(solution, X, Y);
[sxxc, sxyc, syxc, syyc, sxx, sxy, syx, syy] = evaluate_stress(solution, X, Y);

% Stresses from velocity gradient (correct)
%velsxx = 2*Uxc;
%velsxy = Vxc+Uyc;
%velsyx = Uyc+Vxc;
%velsyy = 2*Vyc;

% convert Cartesian velocity to radial and angular velocity, using 
% relationships:
% theta_hat = -i sin(theta) +j cos(theta) = i y/r + j x/r
% r_hat = i cos(theta) + j sin(theta) = i x/r + j y/r
Ur = (Uc.*X + Vc.*Y)./sqrt(X.^2 + Y.^2);
Utheta = (-Uc.*Y + Vc.*X)./sqrt(X.^2 + Y.^2);

% quantities for exact solution
ro = input_params.radii(1);
ri = input_params.radii(2);
omega = input_params.omega;
A = -omega/(ro^2*(1/ri^2-1/ro^2));
B = omega/(1/ri^2-1/ro^2);

% exact solution in polar coordinates
exact_solution_r = @(x,y) zeros(size(x));
exact_solution_theta = @(x,y) A*sqrt(x.^2 + y.^2) + B./sqrt(x.^2 + y.^2);
exact_solution_ux = @(x,y) (2*B*x.*y)./(x.^2+y.^2).^2;
exact_solution_uy = @(x,y) B*(y.^2-x.^2)./(x.^2+y.^2).^2 - A;
exact_solution_vx = @(x,y) B*(y.^2-x.^2)./(x.^2+y.^2).^2 + A;
exact_solution_vy = @(x,y) -(2*B*x.*y)./(x.^2+y.^2).^2;
exact_solution_avg_velocity = pi*(A*(ro^2-ri^2)+2*B*log(ro/ri));
exact_solution_vorticity = @(x,y) 2*A*ones(length(x),length(y));
sigma = exact_solution_sigma(X,Y,A,B);

%% on-surface
[usurf,vsurf] = evaluate_velocity_on_surface(solution,solution);
x = real(problem.domain.z);
y = imag(problem.domain.z);
nu1 = real(problem.domain.zp)./abs(problem.domain.zp);
nu2 = imag(problem.domain.zp)./abs(problem.domain.zp);
n1 = real(-1i*problem.domain.zp)./abs(problem.domain.zp);
n2 = imag(-1i*problem.domain.zp)./abs(problem.domain.zp);
ursurf = usurf.*n1 + vsurf.*n2;
uthetasurf = usurf.*nu1 + vsurf.*nu2;

%%
if test
    h = figure();
    
    subplot(5,2,1);
    plot_domain(problem, h);
    hold on;
    contourf(X,Y, log10(abs(Ur - exact_solution_r(X,Y))+eps));
    colorbar;
    axis equal
    title('u_r: log_{10}(error in radial velocity)');
    
    subplot(5,2,2);    
    plot_domain(problem, h);
    hold on;
    contourf(X,Y,log10(abs(Utheta - exact_solution_theta(X,Y))./...
        max(max(abs(exact_solution_theta(X,Y)))) + eps));
    colorbar;
    axis equal
    title('u_{phi}: log_{10}(relative error angular velocity)');
    
    subplot(5,2,3);
    plot_domain(problem, h);
    hold on
    contourf(X,Y,Pc);
    colorbar
    axis equal
    title('pressure');
    
    subplot(5,2,4);
    plot_domain(problem, h);
    hold on
    contourf(X,Y,log10(abs(omegac - exact_solution_vorticity(X,Y))./...
        max(max(abs(exact_solution_vorticity(X,Y)))) + eps));
    colorbar
    axis equal
    title('vorticity: log_{10}(relative error)');
    
    subplot(5,2,5);
    plot_domain(problem, h);
    hold on
    contourf(X,Y,log10(abs(Pxc) + eps));
    colorbar
    axis equal
    title('dp/dx: log_{10}(absolute error)');

    subplot(5,2,6);
    plot_domain(problem, h);
    hold on
    contourf(X,Y,log10(abs(Pyc) + eps));
    colorbar
    axis equal
    title('dp/dy: log_{10}(absolute error)');

    subplot(5,2,7);
    plot_domain(problem, h);
    hold on
    contourf(X,Y,log10(abs(Uxc - exact_solution_ux(X,Y))./...
        max(max(abs(exact_solution_ux(X,Y)))) + eps));
    colorbar
    axis equal
    title('du/dx: log_{10}(relative error)');
    
    subplot(5,2,8);
    plot_domain(problem, h);
    hold on
    contourf(X,Y,log10(abs(Uyc - exact_solution_uy(X,Y))./...
        max(max(abs(exact_solution_uy(X,Y)))) + eps));
    colorbar
    axis equal
    title('du/dy: log_{10}(relative error)');

    subplot(5,2,9);
    plot_domain(problem, h);
    hold on
    contourf(X,Y,log10(abs(Vxc - exact_solution_vx(X,Y))./...
        max(max(abs(exact_solution_vx(X,Y)))) + eps));
    colorbar
    axis equal
    title('dv/dx: log_{10}(relative error)');
    
    subplot(5,2,10);
    plot_domain(problem, h);
    hold on
    contourf(X,Y,log10(abs(Vyc - exact_solution_vy(X,Y))./...
        max(max(abs(exact_solution_vy(X,Y)))) + eps));
    colorbar
    axis equal
    title('dv/dy: log_{10}(relative error)');
    
    h2 = figure();
    subplot(2,2,1);
    plot_domain(problem, h2);
    hold on;
    contourf(X,Y,log10(abs(sxxc - sigma(:,:,1))./...
        max(max(abs(sigma(:,:,1)))) + eps));
    colorbar;
    axis equal;
    title('sxx: log_{10} relative error');

    subplot(2,2,2);
    plot_domain(problem, h2);
    hold on;
    contourf(X,Y,log10(abs(sxyc - sigma(:,:,3))./...
        max(max(abs(sigma(:,:,3)))) + eps));
    colorbar;
    axis equal;
    title('sxy: log_{10} relative error');

    subplot(2,2,3);
    plot_domain(problem, h2);
    hold on;
    contourf(X,Y,log10(abs(syxc - sigma(:,:,2))./...
        max(max(abs(sigma(:,:,2)))) + eps));
    colorbar;
    axis equal;
    title('syx: log_{10} relative error');

    subplot(2,2,4);
    plot_domain(problem, h2);
    hold on;
    contourf(X,Y,log10(abs(syyc - sigma(:,:,4))./...
        max(max(abs(sigma(:,:,4)))) + eps));
    colorbar;
    axis equal;
    title('syy: log_{10} relative error');
end

function S = exact_solution_sigma(X,Y,A,B)
[m,n] = size(X);
S = zeros(m,n,4);
for i = 1:m
    for j = 1:n
        x = X(i,j);
        y = Y(i,j);

        exact_solution_ux = (2*B*x.*y)./(x.^2+y.^2).^2;
        exact_solution_uy = B*(y.^2-x.^2)./(x.^2+y.^2).^2 - A;
        exact_solution_vx = B*(y.^2-x.^2)./(x.^2+y.^2).^2 + A;
        exact_solution_vy = -(2*B*x.*y)./(x.^2+y.^2).^2;

        U = [exact_solution_ux exact_solution_vx;
             exact_solution_uy exact_solution_vy];
        s = U + U';
        S(i,j,:) = s(:);
    end
end
end
