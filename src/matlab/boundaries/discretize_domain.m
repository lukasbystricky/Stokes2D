function domain = discretize_domain(geometries, panels, varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Discretizes the solid boundaries.
%
% Input:
%       geometries, functions describing the shape of the walls, cell array
%           of length 1x number of walls
%       panels, number on each wall, vector of length 1x number of walls
%       Lx, periodic length in x direction
%       Ly, periodic length in y direction
%
% Output:
%       domain, struct containt solid G-L discretisation points, 
%       derivatives, G-L weights, and gridding information for special 
%       quadrature
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% if domain is periodic user should provide Lx, Ly
periodic = 0;
if nargin == 4
    periodic = 1;
    Lx = varargin{1};
    Ly = varargin{2};
end

% the user may wish to provide the centers of solid, for completion flow
if nargin == 3
    centers = varargin{1};
else
    centers = [];
end

%% Allocate arrays
% discretized values using 16 points per panel
z = zeros(sum(panels)*16,1);     
zp = zeros(sum(panels)*16,1);     
zpp = zeros(sum(panels)*16,1);     
quad_weights = zeros(sum(panels)*16,1);      

% upsampled 32 point discretization
z32 = zeros(sum(panels)*32,1);     
zp32 = zeros(sum(panels)*32,1);   
zpp32 = zeros(sum(panels)*32,1); 
quad_weights32 = zeros(sum(panels)*32,1);     

% indices of points on wall
wall_indices = zeros(length(geometries),2);

% mapping of points to wall number
pts2wall = zeros((sum(panels))*16,1); 

% mapping of panels to wall number
panels2wall = zeros(sum(panels),1); 

% panel breakpoints in physical space
panel_breaks = zeros(sum(panels) + length(geometries),1);

% endpoints of each panel in physical space
panel_endpoints = zeros(sum(panels),2);

%% Discretize each wall using Gauss-Legendre quadrature
ptr = 1;
ptr32 = 1;
ptrPan = 1;
ptrPanBreak = 1;

for j = 1:length(geometries)
   
    % get 16 and 32 point quadrature nodes/weights for wall j
    % the walls are discretized from 0 to 2*pi in parameter space
    [Tn,Wn] = GLinterval(0,2*pi,panels(j));
    [Tn32,Wn32] = GLinterval32(0,2*pi,panels(j));
    
    % discretize wall j
    [zn,zpn,zppn] = geometries{j}(Tn);
    [zn32,zpn32,zppn32] = geometries{j}(Tn32);
    
    % add to arrays
    z(ptr:ptr+16*panels(j)-1) = zn;
    zp(ptr:ptr+16*panels(j)-1) = zpn;
    zpp(ptr:ptr+16*panels(j)-1) = zppn;
    z32(ptr32:ptr32+32*panels(j)-1) = zn32;
    zp32(ptr32:ptr32+32*panels(j)-1) = zpn32;
    zpp32(ptr32:ptr32+32*panels(j)-1) = zppn32;
    
    quad_weights(ptr:ptr+16*panels(j)-1) = Wn;
    quad_weights32(ptr32:ptr32+32*panels(j)-1) = Wn32;
    
    % mapping of points to wall
    pts2wall((ptr-1)+(1:16*panels(j))) = j-1; 
    panels2wall(ptrPan:ptrPan+panels(j)-1) = ones(panels(j),1)*(j-1);
    
    % edges of panels in physical space
    panel_breaks_tmp =  geometries{j}(linspace(0,2*pi,panels(j)+1)');
    panel_breaks(ptrPanBreak: ptrPanBreak + panels(j)) = ...
            panel_breaks_tmp;
    
    panel_endpoints(ptrPan:ptrPan + panels(j)-1,1) = panel_breaks_tmp(1:end-1);
    panel_endpoints(ptrPan:ptrPan + panels(j)-1,2) = panel_breaks_tmp(2:end);
    
    % panel endpoints of wall j
    wall_indices(j,:) = [ptr ptr+16*panels(j)-1];
   
    ptr = ptr + 16*panels(j);
    ptr32 = ptr32 + 32*panels(j);
    ptrPan = ptrPan + panels(j);
    ptrPanBreak = ptrPanBreak + panels(j) + 1;
end

panel_lengths = abs(panel_endpoints(:,2)-panel_endpoints(:,1));
mean_panel_length = mean(panel_lengths);

%% Grid domain, this is useful for special quadrature

% for periodic domains, first box points in reference cell plus points in 
% neighbouring cells

if periodic
    xmin = -1.5*Lx;
    ymin = -1.5*Ly;
    
    XBoxes = 3*ceil(Lx/mean_panel_length);
    YBoxes = 3*ceil(Ly/mean_panel_length);
else
    xmin = min(real(z))- mean_panel_length/2;
    ymin = min(imag(z))- mean_panel_length/2;
    
    xmax = max(real(z)) + mean_panel_length/2;
    ymax = max(real(z)) + mean_panel_length/2;
    
    Lx = xmax - xmin;
    Ly = ymax - ymin;
    
    XBoxes = ceil(Lx/mean_panel_length);
    YBoxes = ceil(Ly/mean_panel_length);
end

gridTmp = cell(XBoxes,YBoxes);

for n = 1:sum(panels)
    z1 = panel_endpoints(n,2);
    z2 = panel_endpoints(n,1);
    mid = (z1+z2)/2;
    zrel = (mid -(xmin+1i*ymin))/mean_panel_length;
    midx = floor(real(zrel))+1;
    midy = floor(imag(zrel))+1;
    radius = ceil(2.5*(panel_lengths(n)/mean_panel_length) + 10*eps);
    for x = midx-radius:midx+radius
        for y = midy-radius:midy+radius
            if x > 0 && x <= XBoxes && y > 0 && y <= YBoxes
                gridTmp{x,y} = [gridTmp{x,y} n];
            end
        end
    end
end

if periodic
    % map everything to reference cell
    domain_grid = cell(XBoxes/3, YBoxes/3);
    
    for i = 1:XBoxes/3
        for j = 1:YBoxes/3
            domain_grid{i,j} = unique([gridTmp{i,j}, ...
                gridTmp{i+XBoxes/3,j}, gridTmp{i+2*XBoxes/3,j},...
                gridTmp{i,j+YBoxes/3}, gridTmp{i,j+2*YBoxes/3},...
                gridTmp{i+XBoxes/3, j+YBoxes/3}, ...
                gridTmp{i+2*XBoxes/3, j+YBoxes/3},...
                gridTmp{i+XBoxes/3, j+2*YBoxes/3},...
                gridTmp{i+2*XBoxes/3, j+2*YBoxes/3}]);
        end
    end
    
    reference_cell = [-Lx/2, Lx/2, -Ly/2, Ly/2];
    
else
    domain_grid = gridTmp;
    
    reference_cell = [xmin, xmax, ymin, ymax];
end

%Save the Solid information in a structure
domain = struct('z',z,'zp',zp,'zpp',zpp,'quad_weights',quad_weights,...
    'wazp',quad_weights.*abs(zp),'z32',z32,'zp32',zp32,'zpp32',zpp32,...
    'quad_weights32',quad_weights32,'wazp32',quad_weights32.*abs(zp32),...
    'wall_indices',wall_indices,'panel_breaks',panel_breaks,...
    'panel_endpoints', panel_endpoints, 'grid', {domain_grid},...
    'reference_cell', reference_cell,...
    'mean_panel_length',mean_panel_length,'pts2wall',pts2wall,...
    'Lx', Lx, 'Ly', Ly, 'centers', centers);

[Nrows,Ncols] = size(domain.grid);
domain.extra.Nrows = Nrows;
domain.extra.Ncols = Ncols;
tmp=reshape(domain.grid,1,Ncols*Nrows);
gridSolidmat = -1*ones(2*length(panel_breaks),Ncols*Nrows);
maxNpans=0;
for j=1:Nrows*Ncols
    Npans=length(tmp{j});
    if Npans>maxNpans
        maxNpans=Npans;
    end
    gridSolidmat(1:Npans,j)=tmp{j}';
end
domain.extra.gridSolidmat=gridSolidmat(1:maxNpans,:);
domain.extra.panels2wall = panels2wall;

% load special quadrature correction for on-surface evaluation of Stokeslet
load 'LmodMatSame.mat'
domain.Lmod = [Lmod;rot90(Lmod,2)];

