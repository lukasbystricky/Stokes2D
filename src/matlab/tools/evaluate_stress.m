function [sxxc, sxyc, syxc, syyc, sxx, sxy, syx, syy, X, Y] = ...
		evaluate_stress(solution, varargin)
%EVALUTATE_STRESS evaluates the stress at either 
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
% -sxxc, sxyc, syxc, syyc: components of the stress evaluated at X,Y
% -sxx, sxy, syx, syy: component of stress without special quadrature 
%		correction
% -X,Y: x and y coordinates of evaluation points

%disp('Evaluating stress...');

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

if solution.problem.periodic

     % works with this code, and with mex_new
    [sigmaxx_slp, sigmayx_slp] = StokesSLP_stress_ewald_2p(xsrc, ysrc, X(:), Y(:),...
        solution.q(:,1).*weights, solution.q(:,2).*weights, ones(ntar,1),...
        zeros(ntar,1), Lx, Ly, 'verbose', 0);

    [sigmaxy_slp, sigmayy_slp] = StokesSLP_stress_ewald_2p(xsrc, ysrc, X(:), Y(:),...
        solution.q(:,1).*weights, solution.q(:,2).*weights, zeros(ntar,1),...
        ones(ntar,1), Lx, Ly, 'verbose', 0);
    
    % using mex_new
%     nsrc = length(xsrc);
%     [sigmaxx_slp, sigmayx_slp] = StokesSLP_stress_ewald_2p(xsrc, ysrc, X(:), Y(:),...
%         solution.q(:,1).*weights, solution.q(:,2).*weights, ones(nsrc,1),...
%         zeros(nsrc,1), Lx, Ly, 'verbose', 0);
% 
%     [sigmaxy_slp, sigmayy_slp] = StokesSLP_stress_ewald_2p(xsrc, ysrc, X(:), Y(:),...
%         solution.q(:,1).*weights, solution.q(:,2).*weights, zeros(nsrc,1),...
%         ones(nsrc,1), Lx, Ly, 'verbose', 0);

    sigma1_slp = sigmaxx_slp + 1i*sigmayx_slp;
    sigma2_slp = sigmaxy_slp + 1i*sigmayy_slp;

    if isinf(solution.problem.eta)
        sigma1 = sigma1_slp;
        sigma2 = sigma2_slp;
    else
        [sigmaxx_dlp, sigmayx_dlp] = StokesDLP_stress_ewald_2p(xsrc, ysrc, X(:), Y(:),...
            n1, n2,solution.q(:,1).*weights, solution.q(:,2).*weights, ones(ntar,1),...
            zeros(ntar,1), Lx, Ly, 'verbose', 0);

        [sigmaxy_dlp, sigmayy_dlp] = StokesDLP_stress_ewald_2p(xsrc, ysrc, X(:), Y(:),...
            n1, n2, solution.q(:,1).*weights, solution.q(:,2).*weights, zeros(ntar,1),...
            ones(ntar,1), Lx, Ly, 'verbose', 0);

        sigma1_dlp = sigmaxx_dlp + 1i*sigmayx_dlp;
        sigma2_dlp = sigmaxy_dlp + 1i*sigmayy_dlp;
        
        sigma1 = solution.problem.eta*sigma1_slp + sigma1_dlp;
        sigma2 = solution.problem.eta*sigma2_slp + sigma2_dlp;
    end
    
    % correct using special quadrature
    % for special quadrature points must be inside reference cell
    Xtar_sq = mod(X+Lx/2,Lx)-Lx/2;
    Ytar_sq = mod(Y+Ly/2,Ly)-Ly/2 + 1e-60;

    % correct using special quadrature
    [sigma1_slp_corrected,~] = mex_SQ_slp_stress(Xtar_sq(:) + 1i*Ytar_sq(:),...
        domain.z, domain.zp, domain.quad_weights, domain.panel_breaks,...
        domain.wazp, domain.z32, domain.zp32, domain.quad_weights32, domain.wazp32, ...
        solution.q(:,1)+1i*solution.q(:,2),ones(ntar,1)+1e-60*1i*ones(ntar,1),...
        sigma1_slp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
        domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
        domain.reference_cell, true);

    [sigma2_slp_corrected,~] = mex_SQ_slp_stress(Xtar_sq(:) + 1i*Ytar_sq(:),...
        domain.z, domain.zp, domain.quad_weights, domain.panel_breaks,...
        domain.wazp, domain.z32, domain.zp32, domain.quad_weights32, domain.wazp32, ...
        solution.q(:,1)+1i*solution.q(:,2),1e-60*ones(ntar,1)+1i*ones(ntar,1),...
        sigma2_slp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
        domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
        domain.reference_cell, true);

    if isinf(solution.problem.eta)
        sigma1c = sigma1_slp_corrected;
        sigma2c = sigma2_slp_corrected;
    else
        % correct using special quadrature
        [sigma1_dlp_corrected,~] = mex_SQ_dlp_stress(Xtar_sq(:) + 1i*Ytar_sq(:),...
            domain.z, domain.zp, domain.quad_weights, domain.panel_breaks,...
            domain.wazp, domain.z32, domain.zp32, domain.quad_weights32, domain.wazp32, ...
            solution.q(:,1)+1i*solution.q(:,2),ones(ntar,1)+1e-60*1i*ones(ntar,1),...
            sigma1_dlp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
            domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
            domain.reference_cell, true);

        [sigma2_dlp_corrected,~] = mex_SQ_dlp_stress(Xtar_sq(:) + 1i*Ytar_sq(:),...
            domain.z, domain.zp, domain.quad_weights, domain.panel_breaks,...
            domain.wazp, domain.z32, domain.zp32, domain.quad_weights32, domain.wazp32, ...
            solution.q(:,1)+1i*solution.q(:,2),1e-60*ones(ntar,1)+1i*ones(ntar,1),...
            sigma2_dlp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
            domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
            domain.reference_cell, true);
        
        sigma1c = solution.problem.eta*sigma1_slp_corrected + sigma1_dlp_corrected;
        sigma2c = solution.problem.eta*sigma2_slp_corrected + sigma2_dlp_corrected;
    end
else % completion formulation
    sigma1_dlp = evaluate_double_layer_stress_direct(X(:),Y(:),...
        domain.z,ones(ntar,1)+1e-60*1i*ones(ntar,1),domain.zp,solution.q(:,1)+1i*solution.q(:,2),domain.quad_weights);
    sigma2_dlp = evaluate_double_layer_stress_direct(X(:),Y(:),...
        domain.z,1e-60*ones(ntar,1)+1i*ones(ntar,1),domain.zp,solution.q(:,1)+1i*solution.q(:,2),domain.quad_weights);
    
    % correct using special quadrature
    [sigma1_dlp_corrected,~] = mex_SQ_dlp_stress(X(:)+1i*(Y(:)+1e-60),...
        domain.z, domain.zp, domain.quad_weights, domain.panel_breaks,...
        domain.wazp, domain.z32, domain.zp32, domain.quad_weights32, domain.wazp32, ...
        solution.q(:,1)+1i*solution.q(:,2),ones(ntar,1)+1e-60*1i*ones(ntar,1),...
        sigma1_dlp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
        domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
        domain.reference_cell, false);

    [sigma2_dlp_corrected,~] = mex_SQ_dlp_stress(X(:)+1i*(Y(:)+1e-60),...
        domain.z, domain.zp, domain.quad_weights, domain.panel_breaks,...
        domain.wazp, domain.z32, domain.zp32, domain.quad_weights32, domain.wazp32, ...
        solution.q(:,1)+1i*solution.q(:,2),1e-60*ones(ntar,1)+1i*ones(ntar,1),...
        sigma2_dlp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
        domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
        domain.reference_cell, false);
    
    % add on completion flow gradients from rotlets and Stokeslets
    [sigma1_S, sigma1_R] = completion_contribution_stress(domain.centers(2:end),...
        X(:)+1i*Y(:),ones(ntar,1)+1i*zeros(ntar,1),solution.forces,solution.torques);
    [sigma2_S, sigma2_R] = completion_contribution_stress(domain.centers(2:end),...
        X(:)+1i*Y(:),zeros(ntar,1)+1i*ones(ntar,1),solution.forces,solution.torques);
    
    sigma1 = sigma1_dlp + sigma1_S + sigma1_R;
    sigma2 = sigma2_dlp + sigma2_S + sigma2_R;
    
    sigma1c = sigma1_dlp_corrected + sigma1_S + sigma1_R;
    sigma2c = sigma2_dlp_corrected + sigma2_S + sigma2_R;
end

if solution.problem.periodic && solution.trim
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

    sigma1(outside) = nan+1i*nan;
    sigma2(outside) = nan+1i*nan;
    sigma1c(outside) = nan+1i*nan;
    sigma2c(outside) = nan+1i*nan;

    X(outside) = nan;
    Y(outside) = nan;
end

sxx = reshape(real(sigma1), size(X));
sxy = reshape(real(sigma2), size(X));
syx = reshape(imag(sigma1), size(X));
syy = reshape(imag(sigma2), size(X));
sxxc = reshape(real(sigma1c), size(X));
sxyc = reshape(real(sigma2c), size(X));
syxc = reshape(imag(sigma1c), size(X));
syyc = reshape(imag(sigma2c), size(X));

% add on driving pressure
if solution.problem.periodic
    sxx = sxx - solution.problem.pressure_gradient_x*X;
    sxx = sxx - solution.problem.pressure_gradient_y*Y;
    syy = syy - solution.problem.pressure_gradient_x*X;
    syy = syy - solution.problem.pressure_gradient_y*Y;
    sxxc = sxxc - solution.problem.pressure_gradient_x*X;
    sxxc = sxxc - solution.problem.pressure_gradient_y*Y;
    syyc = syyc - solution.problem.pressure_gradient_x*X;
    syyc = syyc - solution.problem.pressure_gradient_y*Y;
end
