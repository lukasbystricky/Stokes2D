close all;
clearvars;
clc;
format long;
rng(42);

% create input structure
input_params = default_input_params('channel_blobs', 1);

% modify structure as needed, or add additional problem-dependent params
Lx = 2;
Ly = 2;
h = 0.5;
npan = 10;

input_params.box_size = [Lx,Ly];
input_params.panels = npan;
input_params.h = h;
input_params.eta = inf;
input_params.plot_domain = 0;
input_params.alpha = -1e-1;
input_params.slip = 1;

% set up radii and centers, centers given as complex numbers (x+iy)
M = 3;              % number of blobs
delta = 0.1;        % minimum distance between blobs/blobs and blobs/walls
rmin = 0.1;         % minimum radius
rmax = 0.2;         % maximium radius

centers = zeros(M,1);
radii = zeros(M,1);

for i = 1:M
    overlap = true;
    
    while overlap
        r = rmin + (rmax-rmin).*rand(1,1);

        xmin = -Lx/2 + r;
        xmax = Lx/2 - r;
        ymin = -h/2 + r;
        ymax = h/2 - r;

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

% initialize problem
problem = channel_blobs_periodic(input_params);
domain = problem.domain;

% normalized tangent
nu1 = real(domain.zp)./abs(domain.zp);
nu2 = imag(domain.zp)./abs(domain.zp);

% normalized normal
n1 = real(-1i*domain.zp)./abs(domain.zp);
n2 = imag(-1i*domain.zp)./abs(domain.zp);

%% solve the problem
N = length(problem.domain.z);
rhs = [zeros(N/2,1); zeros(N/2,1); zeros(N/2,1); zeros(N/2,1); 
       problem.pressure_gradient_x; problem.pressure_gradient_y];
problem.alpha = 1e-1*[ones(M*npan*16,1); -ones(16*npan,1); -ones(16*npan,1)];
solution = solve_stokes(problem,rhs,'fmm',0);
solution.local_indices = 1:length(solution.q);

%%
[U,V,X,Y] = evaluate_velocity(solution, 300, 'fmm', 0, 'verbose', 0);
[Uq,Vq,Xq,Yq] = evaluate_velocity(solution, 50, 'fmm', 0, 'verbose', 0);
[Usurf,Vsurf] = evaluate_velocity_on_surface(solution, solution);

%%
% normal velocity (u \dot n)
udotn = Usurf.*n1 + Vsurf.*n2;

% tangential velocity (u \dot \nu)
udotnu = Usurf.*nu1 + Vsurf.*nu2;

%%
figure;
subplot(1,2,1);
contourf(X,Y,U,20);
colorbar;
axis equal;
title('u');

subplot(1,2,2);
contourf(X,Y,V,20);
colorbar;
axis equal;
title('v');

figure;
quiver(Xq,Yq,Uq,Vq,2);
hold on;
for i = 1:M+2
    plot(domain.z(domain.wall_indices(i,1):domain.wall_indices(i,2)),'k')
end
axis equal;
xlim([-Lx/2,Lx/2]);
ylim([-h,h]);
title('velocity field');

figure;
subplot(1,2,1);
plot(Usurf);
grid on;
title('u surf');

subplot(1,2,2);
plot(Vsurf);
grid on;
title('v surf');

figure;
subplot(1,2,1);
plot(udotn);
grid on;
title('u dot n');

subplot(1,2,2);
plot(udotnu);
grid on;
title('u dot nu');
