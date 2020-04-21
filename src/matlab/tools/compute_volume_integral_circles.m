function [I1, I2, X, Y, U1, U2] = compute_volume_integral_circles...
                    (centers, radii, solution, varargin)
% COMPUTE_VOLUME_INTEGRALS_CIRCLES computes the volume integrals of the x
% and y components of a velocity field over a domain that is a box with
% circlular holes cut out. 
%
% inputs:
% - centers: vector of complex numbers denoting the centers of each circle
% - radii: vector of radii
% - solution: solution structure, containing domain information and BIE
% solution
% - varargin: optional parameters, specified as a name-value pair. Can be
% any of:
%       - hmax: maximum mesh spacing of volume grid (default 0.3)
%       - N: related to the number of quadrature points in each volume
%       element, the number of points is N^2



%% read in optional input parameters

quad_pts = 3; % default values
hmax = 0.3;

if nargin > 3
    % Go through all other input arguments and assign parameters
    jv = 1;
    while jv <= length(varargin)-1
       switch varargin{jv}
  
           case 'hmax'
               hmax = varargin{jv+1};
               
           case 'N'
               quad_pts = varargin{jv+1};
       end
       jv = jv + 2;
    end
end

%% create domain

Lx = solution.problem.domain.Lx;
Ly = solution.problem.domain.Ly;
box = [3; 4; -Lx/2; -Lx/2; Lx/2; Lx/2; -Ly/2; Ly/2; Ly/2; -Ly/2];

n_circles = length(radii);
geo_names = cell(9*n_circles+1,1);
geo_names{1} = 'box';

circles = zeros(10, 9*n_circles);

% create circles, and periodic replicates
for i = 1:n_circles
    
    for j = -1:1
        for k = -1:1
            circles(1,9*(i-1)+3*(j+1)+k+2) = 1;
            circles(2,9*(i-1)+3*(j+1)+k+2) = real(centers(i)) + j*Lx;
            circles(3,9*(i-1)+3*(j+1)+k+2) = imag(centers(i)) + k*Lx;
            circles(4,9*(i-1)+3*(j+1)+k+2) = radii(i);
            
            geo_names{9*(i-1)+3*(j+1)+k+3} = ...
                        ['circle',num2str(9*(i-1)+3*(j+1)+k+2)];
        end
    end
    
end

gd = [box, circles];
ns = char(geo_names)';

sf = 'box';

for i = 1:9*n_circles
    sf = strcat(sf, ['-circle', num2str(i)]);
end

[dl,~] = decsg(gd,sf,ns);

pdegplot(dl,'EdgeLabels','on','FaceLabels','on')
axis equal

%% create mesh
model = createpde;
geometryFromEdges(model,dl);

mesh = generateMesh(model, 'Hmax', hmax, 'geometricOrder', 'linear');

hold on
pdemesh(model)
drawnow

%% collect quadrature points and weights

n_elements = size(mesh.Elements,2);

quad_box = zeros(n_elements * quad_pts^2, 2);
w_quad_box = zeros(n_elements * quad_pts^2, 1);

for i = 1:n_elements
    tri_indices = mesh.Elements(:,i);
    v = mesh.Nodes(:,tri_indices);
    
    [X,Y,Wx,Wy]=triquad(quad_pts,v');
    
    quad_box((i-1)*quad_pts^2 + 1:i*quad_pts^2, 1) = X(:);
    quad_box((i-1)*quad_pts^2 + 1:i*quad_pts^2, 2) = Y(:);
    
    tmp = Wx*Wy';
    
    w_quad_box((i-1)*quad_pts^2 + 1:i*quad_pts^2) = tmp(:);
end

%% compute integral
X = quad_box(:,1);
Y = quad_box(:,2);
[U1, U2] = evaluate_velocity(solution, X, Y);

U1(isnan(U1)) = 0;
U2(isnan(U2)) = 0;

I1 = sum(U1.*w_quad_box);
I2 = sum(U2.*w_quad_box);


