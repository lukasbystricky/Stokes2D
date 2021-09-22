% Computes the velocity between two rotating cylinders, i.e. Taylor-Couette
% flow. The outer cylinder is fixed, while the inner cylinder rotates
% counter clockwise with angular velocity omega. Then, the exact solution
% is known, with angular velocity given by A*r + B/r, where
% A = omega r_out^2 / (r_out^2 - r_in^2),
% B = -omega r_out^2 r_in^2 /(r_out^2 - r_in^2)

close all
clearvars
clc

test = 1;

% create input structure
input_params = default_input_params('couette', 0);

% modify structure as needed, or add additional problem-dependent params
input_params.panels = [20, 10];
input_params.plot_domain = 0;

% set up radii and centers, centers given as complex numbers (x+iy)
input_params.radii = [2, 0.5]; % the outer wall must be the first radii
input_params.centers = [0, 0];
input_params.omega = 1; % angular velocity of outer wall

problem = couette(input_params);

% solve the problem
solution = solve_stokes(problem);

% grid in polar coordinates with M number of radial and angular points
M = 200;
r = linspace(input_params.radii(2)+1e-1, input_params.radii(1) - 1e-1, M);
h = (2*pi)/M;
theta = (0:h:2*pi);
[R, T] = meshgrid(r, theta);

X = R.*cos(T);
Y = R.*sin(T);

[Uc, Vc, X, Y, U, V] = evaluate_velocity(solution, X, Y);
Pc = evaluate_pressure(solution, X, Y);

% convert Cartesian velocity to radial and angular velocity, using 
% relationships:
% theta_hat = -i sin(theta) +j cos(theta) - i y/r + j x/r
% r_hat = i cos(theta) + j sin(theta) = i x/r + j y/r
Ur = (Uc.*X + Vc.*Y)./sqrt(X.^2 + Y.^2);
Utheta = (-Uc.*Y + Vc.*X)./sqrt(X.^2 + Y.^2);

% Quantities for exact solution
ro = input_params.radii(1);
ri = input_params.radii(2);
omega = input_params.omega;
A = -omega/(ro^2*(1/ri^2-1/ro^2));
B = omega/(1/ri^2-1/ro^2);

 % exact solution in polar coordinates
exact_solution_r = @(x,y) zeros(size(x));
exact_solution_theta = @(x,y) A*sqrt(x.^2 + y.^2) + B./sqrt(x.^2 + y.^2);
exact_solution_pressure = @(x,y) A^2*(x.^2+y.^2)/2 + 2*A*B*log(sqrt(x.^2+y.^2)) - B^2./(x.^2+y.^2);

h = figure();

if test
    
    subplot(2,2,1);
    plot_domain(problem, h);
    hold on;
    contourf(X,Y, log10(abs(Ur - exact_solution_r(X,Y))+eps));
    colorbar;
    axis equal
    title('log_{10}(error in radial velocity)');
    
    subplot(2,2,2);    
    plot_domain(problem, h);
    hold on;
    contourf(X,Y,log10(abs(Utheta - exact_solution_theta(X,Y))./...
        max(max(abs(exact_solution_theta(X,Y)))) + eps));
    colorbar;
    axis equal
    title('log_{10}(relative error angular velocity)');
    
    subplot(2,2,3);
    plot_domain(problem, h);
    hold on
    contourf(X,Y,Pc);
    colorbar
    axis equal
    title('pressure (computed solution)');

    subplot(2,2,4);
    plot_domain(problem, h);
    hold on
    contourf(X,Y,exact_solution_pressure(X,Y));
    colorbar
    axis equal
    title('pressure (analytical solution)');
else
    
    subplot(2,2,1);
    
    plot_domain(problem, h);
    hold on
    contourf(X,Y,Ur);
    colorbar
    axis equal
    title('radial velocity');
    
    subplot(1,2,2);
    
    plot_domain(problem, h);
    hold on
    contourf(X,Y,Utheta);
    colorbar
    axis equal
    title('angular velocity');
end