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
input_params.nbr_neighbor_pts = 4;
input_params.eta = inf;
input_params.plot_domain = 0;
input_params.alpha = 1e-1;
input_params.slip = 1;
input_params.gmres_tol = 1e-12;

% setup domain
input_params.centers = [-0.6+0.15*1i -0.6-0.25*1i 0.1+0.2*1i 0.72-0.25*1i 0.7+0.15*1i 0.1-0.2*1i];
M = length(input_params.centers);

% initialize problem
problem = poster_blobs_periodic(input_params);
problem.alpha = input_params.alpha*[ones(M*npan*16,1); -ones(16*npan,1); -ones(16*npan,1)];
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

solution = solve_stokes(problem,rhs,'fmm',0);
solution.local_indices = 1:length(solution.q);

%%
[U,V,X,Y] = evaluate_velocity(solution, 400, 'fmm', 0, 'verbose', 0);

%%
R = [29 91 170 214 212 161]./255*1.15;
G = [51 118 179 176 123 63]./255*1.11;
B = [101 161 194 155 93 54]./255*1.15;

val = linspace(0,1,100);

close all;
figure;
plot(problem.domain.z,'k.');
hold on;
surfc(X,Y,abs(U+1i*V));
shading interp;
colorbar;
cmap = colormap([R(:), G(:), B(:)]);
hsv = rgb2hsv(cmap);
cm_data=interp1(linspace(0,1,size(cmap,1)),hsv,val);
cm_data=hsv2rgb(cm_data);
colormap(cm_data);
%colormap('bone');
axis off;
axis equal;

