function domain = discretize_domain_dyad_refine(geometries, panels, nsub, varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Discretizes the solid boundaries.
%
% Input:
%       geometries, functions describing the shape of the walls, cell array
%           of length 1x number of walls
%       panels, number on each wall, vector of length 1x number of walls
%       nsub, number of subdivisions of the first and last panel (seen to
%       the parametrization)
%       Lx, periodic length in x direction
%       Ly, periodic length in y direction
%
% Output:
%       domain, struct containt solid G-L discretisation points, 
%       derivatives, G-L weights, and gridding information for special 
%       quadrature
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% the user may wish to provide the centers of solid, for completion flow
if nargin == 5
    centers = varargin{1};
else
    centers = [];
end

%% Allocate arrays
nq = 16;
nbr_panels = sum(panels)+length(panels)*2*nsub;

% discretized values using nq points per panel
z = zeros(nbr_panels*nq,1);   
zp = zeros(nbr_panels*nq,1);     
zpp = zeros(nbr_panels*nq,1);     
quad_weights = zeros(nbr_panels*nq,1);      
t = zeros(nbr_panels*nq,1);

% upsampled 32 point discretization
z32 = zeros(nbr_panels*32,1);   
zp32 = zeros(nbr_panels*32,1);     
zpp32 = zeros(nbr_panels*32,1);     
quad_weights32 = zeros(nbr_panels*32,1);      

% indices of points on wall
wall_indices = zeros(length(geometries),2);

% mapping of points to wall number
pts2wall = zeros((nbr_panels)*nq,1); 

% mapping of panels to wall number
panels2wall = zeros(nbr_panels,1); 

% panel breakpoints in physical space
panel_breaks = zeros(nbr_panels + length(geometries),1);

% endpoints of each panel in physical space
panel_endpoints = zeros(nbr_panels,2);

%% Discretize each wall using Gauss-Legendre quadrature
ptr = 1;
ptr32 = 1;
ptrPan = 1;
ptrPanBreak = 1;

for j = 1:length(geometries)

    % fine paramatrization
    npan = panels(j);
    npanfin = npan+2*nsub;
    pedges = zeros(npanfin+1,1);
    pedges(1:npan+1) = linspace(0,2*pi,npan+1);
    for k = 1:nsub
        pedges(3:end) = pedges(2:end-1);
        pedges(2) = (pedges(1)+pedges(2))/2;   
    end
    pedges(end-nsub:end) = 2*pi-flipud(pedges(1:nsub+1));

    [Tn,Wn] = GLinterval(0,2*pi,npanfin,pedges);
    [Tn32,Wn32] = GLinterval32(0,2*pi,npanfin,pedges);

    % discretize wall j
    [zn,zpn,zppn] = geometries{j}(Tn);
    [zn32,zpn32,zppn32] = geometries{j}(Tn32);
    
    % add to arrays
    z(ptr:ptr+nq*npanfin-1) = zn;
    zp(ptr:ptr+nq*npanfin-1) = zpn;
    zpp(ptr:ptr+nq*npanfin-1) = zppn;
    z32(ptr32:ptr32+32*npanfin-1) = zn32;
    zp32(ptr32:ptr32+32*npanfin-1) = zpn32;
    zpp32(ptr32:ptr32+32*npanfin-1) = zppn32;
    
    quad_weights(ptr:ptr+nq*npanfin-1) = Wn;
    quad_weights32(ptr32:ptr32+32*npanfin-1) = Wn32;
    t(ptr:ptr+nq*npanfin-1) = Tn;
    
    % mapping of points to wall
    pts2wall((ptr-1)+(1:nq*npanfin)) = j-1;
    panels2wall(ptrPan:ptrPan+npanfin-1) = ones(npanfin,1)*(j-1);
    
    % edges of panels in physical space
    panel_breaks_ttmp = linspace(0,2*pi,npanfin+1)';
    panel_breaks_tmp =  geometries{j}(pedges);
    
    panel_breaks_t(ptrPanBreak: ptrPanBreak + npanfin) = ...
            panel_breaks_ttmp;
    panel_breaks(ptrPanBreak: ptrPanBreak + npanfin) = ...
            panel_breaks_tmp;
    
    panel_endpoints(ptrPan:ptrPan + npanfin-1,1) = panel_breaks_tmp(1:end-1);
    panel_endpoints(ptrPan:ptrPan + npanfin-1,2) = panel_breaks_tmp(2:end);
    
    % panel endpoints of wall j
    wall_indices(j,:) = [ptr ptr+nq*npanfin-1];
   
    ptr = ptr + nq*npanfin;
    ptr32 = ptr32 + 32*npanfin;
    ptrPan = ptrPan + npanfin;
    ptrPanBreak = ptrPanBreak + npanfin + 1;
end

panel_lengths = abs(panel_endpoints(:,2)-panel_endpoints(:,1));
mean_panel_length = mean(panel_lengths);

%% Grid domain, this is useful for special quadrature
xmin = min(real(z))- 10*mean_panel_length;
ymin = min(imag(z))- 10*mean_panel_length;

xmax = max(real(z)) + 10*mean_panel_length;
ymax = max(imag(z)) + 10*mean_panel_length;

Lx = xmax - xmin;
Ly = ymax - ymin;

XBoxes = ceil(Lx/mean_panel_length);
YBoxes = ceil(Ly/mean_panel_length);

gridTmp = cell(XBoxes,YBoxes);
for n = 1:nbr_panels
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

domain_grid = gridTmp;
    
reference_cell = [xmin, xmax, ymin, ymax];

domain = struct('z',z,'zp',zp,'zpp',zpp,'quad_weights',quad_weights,...
    'wazp',quad_weights.*abs(zp),'z32',z32,'zp32',zp32,'zpp32',zpp32,...
    'quad_weights32',quad_weights32,'wazp32',quad_weights32.*abs(zp32),...
    'theta', t, 'wall_indices',wall_indices,'panel_breaks',panel_breaks,...
    'panel_endpoints', panel_endpoints, 'panel_breaks_t', panel_breaks_t,...
    'grid', {domain_grid},'reference_cell', reference_cell,...
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
