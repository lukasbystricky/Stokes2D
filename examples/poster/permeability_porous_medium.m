close all;
clearvars;
clc;
format long;

%rng(33);
rng(2);

% create input structure
input_params = default_input_params('channel_blobs', 1);

% modify structure as needed, or add additional problem-dependent params
Lx = 1;
Ly = 1;
h = 0.5;
npan = 20;

input_params.box_size = [Lx,Ly];
input_params.panels = npan;
input_params.h = h;
input_params.nbr_neighbor_pts = 16;
input_params.gmres_tol = 1e-12;
input_params.eta = inf;
input_params.plot_domain = 0;
input_params.alpha = 0;
input_params.slip = 0;

% set up radii and centers, centers given as complex numbers (x+iy)
M = 20;              % number of blobs
M = 2;
delta = 0.19;
rmin = 0.05;
rmax = 0.15;

centers = zeros(M,1);
radii = zeros(M,1);

for i = 1:M
    overlap = true;
    r = rmin + (rmax-rmin).*rand(1,1);
    
    while overlap
        r = 0.99*r;

        xmin = -Lx/2 + r + delta;
        xmax = Lx/2 - r - delta;
        ymin = -Ly/2 + r + delta;
        ymax = Ly/2 - r - delta;

        cx = xmin + (xmax-xmin).*rand(1,1);
        cy = ymin + (ymax-ymin).*rand(1,1);
        c = cx + 1i*cy;

        for j = 1:M
            if M == 1
                overlap = false;
                break;
            end
            if i ~= j
                if abs(centers(j)-c) <= (radii(j)+r+delta)
                    overlap = true;
                    break;
                else
                    overlap = false;
                end
            end
        end
    end
    
    radii(i) = r;
    centers(i) = cx + 1i*cy;
    
end
input_params.radii = radii;
input_params.centers = centers;

% permeability tensor
K = zeros(2,2);

%% solve for K_{11} and K_{12} by imposing pressure gradient in x direction 
input_params.pressure_drop_x = 1;
input_params.pressure_drop_y = 0;

problem = channel_blobs_periodic(input_params);

rhs = [zeros(2*length(problem.domain.z),1);
       problem.pressure_gradient_x;
       problem.pressure_gradient_y];
   
% solve the problem
solution_x = solve_stokes(problem,rhs);
solution_x.local_indices = 1:length(solution_x.q);

% average velocity is computed already!
%%
K(:,1) = solution_x.u_avg;
%u_avg_x = compute_average_velocity(solution_x,solution_x,-problem.Lx/2,problem.Lx/2,-problem.Ly/2,problem.Ly/2,100,100);
%K(:,1) = u_avg_x;

%% solve for K_{21} and K_{22} by imposing pressure gradient in y direction
input_params.pressure_drop_x = 0;
input_params.pressure_drop_y = 1;

problem = channel_blobs_periodic(input_params);

rhs = [zeros(2*length(problem.domain.z),1);
       problem.pressure_gradient_x;
       problem.pressure_gradient_y];
   
% solve the problem
solution_y = solve_stokes(problem,rhs);
solution_y.local_indices = 1:length(solution_y.q);

% average velocity is computed already!

%%
K(:,2) = solution_y.u_avg
%u_avg_y = compute_average_velocity(solution_y,solution_y,-problem.Lx/2,problem.Lx/2,-problem.Ly/2,problem.Ly/2,100,100);
%K(:,2) = u_avg_y

%% Plot
% Off-surface velocity
[U,V,X,Y] = evaluate_velocity(solution_y, 200, 'fmm', 0, 'verbose', 0);

%%
figure;
hold on;
surfc(X,Y,abs(U+1i*V));
shading interp;
colorbar;
axis equal;
