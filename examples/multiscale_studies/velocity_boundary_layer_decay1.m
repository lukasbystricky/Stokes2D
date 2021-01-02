% Investigation of the decay of the velocity boundary layer for flow over a
% porous bed.

close all
clearvars
clc

% create input structure
input_params = default_input_params('decay_study1', 1);

% modify structure as needed, or add additional problem-dependent params
Lx = 10;
Ly = 10;

input_params.box_size = [Lx,Ly];
input_params.panels = 10;
input_params.plot_domain = 0;
input_params.test = 0;
input_params.pressure_drop_x = 0.1;

% set up radii and centers, centers given as complex numbers (x+iy)
n_layers = 5;
L = input_params.box_size(2)/8;
c = 0.2;% concentration
h = L / n_layers; % spacing between circle centers
x_centers = linspace(h/2, L-h/2, n_layers)';
y_centers = linspace(h/2, L-h/2, n_layers)';
[X,Y] = meshgrid(x_centers, y_centers);

radii = L*sqrt(c/pi)/n_layers*ones(n_layers*n_layers, 1);
centers = X(:) + 1i*Y(:);
input_params.radii = radii;
input_params.centers = 1i*centers;

problem = boundary_layer_test(input_params);

% solve the problem
solution = solve_stokes(problem);


if input_params.test
    close all
    [Uc, Vc, X, Y,U] = evaluate_velocity(solution, 100);

    plot_domain(problem);
    hold on
    
    exact_solution = @(x,y) y;
    contourf(X,Y, log10(abs((Uc - exact_solution(X,Y)+eps))));
    colorbar
    axis equal
    title('log_{10}(relative error)');
else
    
    x = linspace(0,Lx)';
    
    % slices above porous layer
    y_above = [0.5;0.6;0.7;0.8;0.9;0.99]*Ly;
    
    % slices below porous layer
    y_below = zeros(n_layers,1);
    y_below(1) = (centers(1)-radii(1))/2;
    for i = 2:n_layers
        y_below(i) = (centers(i) + centers(i-1))/2;
    end
    
    y = [y_below; y_above];
    
    [X, Y] = meshgrid(x,y);
    
    [Uc, Vc, X, Y] = evaluate_velocity(solution, 1000);

    subplot(1,2,1)
    [C,h] = contourf(X,Y,Uc);
    axis equal
    colorbar;
    %clabel(C,h)
    
    subplot(1,2,2)
    [C,h] = contourf(X,Y,Vc);
    axis equal
    colorbar;
    %clabel(C,h)
    
end
