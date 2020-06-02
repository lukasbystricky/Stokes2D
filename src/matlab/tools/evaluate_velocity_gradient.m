function [Ux, Uy, Vx, Vy, X, Y] = evaluate_velocity_gradient(solution, varargin)
%EVALUTATE_VELOCITY evaluates the pressure at either regular grid over the 
%reference cell or at specified target points. Uses the stresslet identity 
%to identify points inside the domain.
%
% inputs:
% -solution: structure containing the following fields
%   -problem: a problem structure containing problem/geometry information
%   -q : density function, vector of size #npts by 2
% variable number of additional arguments, one of
% 1) N: number of points in each direction
% 2) X and Y, vectors (or matrices) of target points

disp('Evaluating velocity gradient...');

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
    x = linspace(min(xsrc), max(xsrc), N);
    y = linspace(min(ysrc), max(ysrc), N);
    
    [X,Y] = meshgrid(x,y);
else % target points are specified
    X = varargin{1};
    Y = varargin{2};
end

ntar = length(X(:));

[ux_slp, vx_slp] = StokesSLP_gradient_ewald_2p(xsrc, ysrc, X(:), Y(:),...
    solution.q(:,1).*weights, solution.q(:,2).*weights, ones(ntar,1),...
    zeros(ntar,1), Lx, Ly, 'verbose', 1);

[uy_slp, vy_slp] = StokesSLP_gradient_ewald_2p(xsrc, ysrc, X(:), Y(:),...
    solution.q(:,1).*weights, solution.q(:,2).*weights, zeros(ntar,1),...
    ones(ntar,1), Lx, Ly, 'verbose', 1);

if isinf(solution.problem.eta)
    ug1 = ux_slp + 1i*uy_slp;
    ug2 = vx_slp + 1i*vy_slp;
else
    [ux_dlp, vx_dlp] = StokesDLP_gradient_ewald_2p(xsrc, ysrc, X(:), Y(:),...
        n1, n2,solution.q(:,1).*weights, solution.q(:,2).*weights, ones(ntar,1),...
        zeros(ntar,1), Lx, Ly, 'verbose', 1);
    
    [uy_dlp, vy_dlp] = StokesDLP_gradient_ewald_2p(xsrc, ysrc, X(:), Y(:),...
        n1, n2, solution.q(:,1).*weights, solution.q(:,2).*weights, zeros(ntar,1),...
        ones(ntar,1), Lx, Ly, 'verbose', 1);
    
    ug1 = solution.problem.eta*(ux_slp + 1i*uy_slp) + ux_dlp + 1i*uy_dlp;
    ug2 = solution.problem.eta*(vx_slp + 1i*vy_slp) + vx_dlp + 1i*vy_dlp;
end

%     [udlp1, udlp2] = StokesDLP_ewald_2p(xsrc, ysrc, X(:), Y(:), n1, n2,...
%         solution.q(:,1).*weights, solution.q(:,2).*weights, Lx, Ly,...
%         'verbose', 1);
%     
%     uslp = uslp1 + 1i*uslp2;
%     udlp = udlp1 + 1i*udlp2;
%     
%     disp('Beginning special quadrature...');
%     
%     % correct using special quadrature
%     [uslp_corrected,~] = mex_SQ_slp(X(:)+1i*Y(:), domain.z, domain.zp,...
%         domain.quad_weights, domain.panel_breaks, domain.wazp, domain.z32,...
%         domain.zp32, domain.quad_weights32, domain.wazp32, ...
%         solution.q(:,1)+1i*solution.q(:,2),...
%         uslp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
%         domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
%         domain.reference_cell);
%     
%     [udlp_corrected,~] = mex_SQ_dlp(X(:)+1i*Y(:), domain.z, domain.zp,...
%         domain.quad_weights, domain.panel_breaks, domain.wazp, domain.z32,...
%         domain.zp32, domain.quad_weights32, domain.wazp32,...
%         solution.q(:,1)+1i*solution.q(:,2),...
%         udlp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
%         domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
%         domain.reference_cell);
%     
%     u_corrected = udlp_corrected + solution.problem.eta*uslp_corrected  + ...
%         solution.u_avg(1) + 1i*solution.u_avg(2);
%     u = udlp + solution.problem.eta*uslp + solution.u_avg(1) +...
%         1i*solution.u_avg(2);
    

% find points inside domain by applying stresslet identity
[test1, test2] = StokesDLP_ewald_2p(xsrc, ysrc, X(:), Y(:), n1, n2,...
        ones(length(n1),1).*weights, zeros(length(n1),1).*weights, Lx, Ly);

    % correct using special quadrature
[test,~] = mex_SQ_dlp(X(:)+1i*Y(:), domain.z, domain.zp, domain.quad_weights, ...
                domain.panel_breaks, domain.wazp, domain.z32, domain.zp32,...
                domain.quad_weights32, domain.wazp32,ones(length(n1),1) + 1e-14*1i,...
                test1 + 1i*test2,domain.mean_panel_length,domain.extra.gridSolidmat, ...
                domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
                domain.reference_cell);

% anything that is greater than 0 is outside the fluid domain
outside = find(solution.problem.stresslet_id_test(real(test)) == 1);
ug1(outside) = nan;
ug2(outside) = nan;

X(outside) = nan;
Y(outside) = nan;

Ux = reshape(real(ug1), size(X));
Uy = reshape(imag(ug1), size(X));
Vx = reshape(real(ug2), size(X));
Vy = reshape(imag(ug2), size(X));