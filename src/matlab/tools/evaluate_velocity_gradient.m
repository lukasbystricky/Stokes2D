function [Uxc, Uyc, Vxc, Vyc, Ux, Uy, Vx, Vy, X, Y] = ...
		evaluate_velocity_gradient(solution, varargin)
%EVALUTATE_VELOCITY_GRADIENT evaluates the velocity gradient  at either 
%regular grid over the reference cell or at specified target points. 
%Uses the stresslet identity to identify points inside the domain.
%
% inputs:
% -solution: structure containing the following fields
%   -problem: a problem structure containing problem/geometry information
%   -q : density function, vector of size #npts by 2
% variable number of additional arguments, one of
% 1) N: number of points in each direction
% 2) X and Y, vectors (or matrices) of target points
% 
% outputs:
% -Uxc, Uyc, Vxc, Vyc: components of the velocity gradient evaluated at X,Y
% -Ux, Uy, Vx, Vy: component of velocity gradient without special quadrature 
%		correction
% -X,Y: x and y coordinates of evaluation points

%disp('Evaluating velocity gradient...');

domain = solution.problem.domain;
if solution.problem.periodic
    Lx = domain.Lx;
    Ly = domain.Ly;
end

xsrc = real(solution.problem.domain.z);
ysrc = imag(solution.problem.domain.z);
n1 = real(-1i*solution.problem.domain.zp)./abs(solution.problem.domain.zp);
n2 = imag(-1i*solution.problem.domain.zp)./abs(solution.problem.domain.zp);
weights = solution.problem.domain.wazp;

if nargin == 2 % N specified, evaluate on regular grid
    N = varargin{1};
    if solution.problem.periodic
        x = linspace(-Lx/2, Lx/2, N);
        y = linspace(-Ly/2, Ly/2, N);
    else
        x = linspace(min(xsrc), max(xsrc), N);
        y = linspace(min(ysrc), max(ysrc), N);
    end
    
    [X,Y] = meshgrid(x,y);
else
    X = varargin{1};
    Y = varargin{2};    
end

ntar = length(X(:));

[ux_slp, vx_slp] = StokesSLP_gradient_ewald_2p(xsrc, ysrc, X(:), Y(:),...
    solution.q(:,1).*weights, solution.q(:,2).*weights, ones(ntar,1),...
    zeros(ntar,1), Lx, Ly, 'verbose', 0);

[uy_slp, vy_slp] = StokesSLP_gradient_ewald_2p(xsrc, ysrc, X(:), Y(:),...
    solution.q(:,1).*weights, solution.q(:,2).*weights, zeros(ntar,1),...
    ones(ntar,1), Lx, Ly, 'verbose', 0);

ug1_slp = ux_slp + 1i*vx_slp;
ug2_slp = uy_slp + 1i*vy_slp;

if isinf(solution.problem.eta)
    ug1 = ug1_slp;
    ug2 = ug2_slp;
else
    [ux_dlp, vx_dlp] = StokesDLP_gradient_ewald_2p(xsrc, ysrc, X(:), Y(:),...
        n1, n2,solution.q(:,1).*weights, solution.q(:,2).*weights, ones(ntar,1),...
        zeros(ntar,1), Lx, Ly, 'verbose', 0);
    
    [uy_dlp, vy_dlp] = StokesDLP_gradient_ewald_2p(xsrc, ysrc, X(:), Y(:),...
        n1, n2, solution.q(:,1).*weights, solution.q(:,2).*weights, zeros(ntar,1),...
        ones(ntar,1), Lx, Ly, 'verbose', 0);
    
    ug1_dlp = ux_dlp + 1i*vx_dlp;
    ug2_dlp = uy_dlp + 1i*vy_dlp;
    
    ug1 = solution.problem.eta*ug1_slp + ug1_dlp;
    ug2 = solution.problem.eta*ug2_slp + ug2_dlp;
end


%disp('Beginning special quadrature...');

% correct using special quadrature
% for special quadrature points must be inside reference cell
Xtar_sq = mod(X+Lx/2,Lx)-Lx/2;
Ytar_sq = mod(Y+Ly/2,Ly)-Ly/2 + 1e-60;

[ug1_slp_corrected,~] = mex_SQ_slp_velocity_grad(Xtar_sq(:) + 1i*Ytar_sq(:), domain.z, domain.zp,...
    domain.quad_weights, domain.panel_breaks, domain.wazp, domain.z32,...
    domain.zp32, domain.quad_weights32, domain.wazp32, ...
    solution.q(:,1)+1i*solution.q(:,2), ones(ntar,1) + 1e-60*1i*ones(ntar,1),...
    ug1_slp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
    domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
    domain.reference_cell,true);

% correct using special quadrature
[ug2_slp_corrected,~] = mex_SQ_slp_velocity_grad(Xtar_sq(:) + 1i*Ytar_sq(:), domain.z, domain.zp,...
    domain.quad_weights, domain.panel_breaks, domain.wazp, domain.z32,...
    domain.zp32, domain.quad_weights32, domain.wazp32, ...
    solution.q(:,1)+1i*solution.q(:,2), 1e-60*ones(ntar,1) + 1i*ones(ntar,1),...
    ug2_slp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
    domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
    domain.reference_cell,true);

if isinf(solution.problem.eta)
    ug1c = ug1_slp_corrected;
    ug2c = ug2_slp_corrected;
else
    
    % correct using special quadrature
    [ug1_dlp_corrected,~] = mex_SQ_dlp_velocity_grad(Xtar_sq(:) + 1i*Ytar_sq(:), domain.z, domain.zp,...
        domain.quad_weights, domain.panel_breaks, domain.wazp, domain.z32,...
        domain.zp32, domain.quad_weights32, domain.wazp32, ...
        solution.q(:,1)+1i*solution.q(:,2), ones(ntar,1) + 1e-60*1i*ones(ntar,1),...
        ug1_dlp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
        domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
        domain.reference_cell,true);
    
    % correct using special quadrature
    [ug2_dlp_corrected,~] = mex_SQ_dlp_velocity_grad(Xtar_sq(:) + 1i*Ytar_sq(:), domain.z, domain.zp,...
        domain.quad_weights, domain.panel_breaks, domain.wazp, domain.z32,...
        domain.zp32, domain.quad_weights32, domain.wazp32, ...
        solution.q(:,1)+1i*solution.q(:,2), 1e-60*ones(ntar,1) + 1i*ones(ntar,1),...
        ug2_dlp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
        domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
        domain.reference_cell,true);
    
    ug1c = solution.problem.eta*ug1_slp_corrected + ug1_dlp_corrected;
    ug2c = solution.problem.eta*ug2_slp_corrected + ug2_dlp_corrected;
end

% find points inside domain by applying stresslet identity
[test1, test2] = StokesDLP_ewald_2p(xsrc, ysrc, X(:), Y(:), n1, n2,...
    ones(length(n1),1).*weights, zeros(length(n1),1).*weights, Lx, Ly);

% correct using special quadrature
[test,~] = mex_SQ_dlp(Xtar_sq(:)+1i*(Ytar_sq(:)+1e-60), domain.z,...
                domain.zp, domain.quad_weights, ...
                domain.panel_breaks, domain.wazp, domain.z32, domain.zp32,...
                domain.quad_weights32, domain.wazp32,ones(length(n1),1) + 1e-14*1i,...
                test1 + 1i*test2,domain.mean_panel_length,domain.extra.gridSolidmat, ...
                domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
                domain.reference_cell,solution.problem.periodic);


% anything that is greater than 0 is outside the fluid domain
outside = find(solution.problem.stresslet_id_test(real(test)) == 1);

ug1(outside) = nan+1i*nan;
ug2(outside) = nan+1i*nan;
ug1c(outside) = nan+1i*nan;
ug2c(outside) = nan+1i*nan;

X(outside) = nan;
Y(outside) = nan;

Ux = reshape(real(ug1), size(X));
Uy = reshape(real(ug2), size(X));
Vx = reshape(imag(ug1), size(X));
Vy = reshape(imag(ug2), size(X));
Uxc = reshape(real(ug1c), size(X));
Uyc = reshape(real(ug2c), size(X));
Vxc = reshape(imag(ug1c), size(X));
Vyc = reshape(imag(ug2c), size(X));
