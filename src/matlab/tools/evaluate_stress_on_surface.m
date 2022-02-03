function [sxx, sxy, syx, syy] = evaluate_stress_on_surface(solution, solution_local, type)
%EVALUTATE_STRESS_ON_SURFACE evaluates the stress at 
%the quadrature points on the surface of a domain. Adds on the jump 
%corresponding to approaching the boundary from the fluid part of the 
%domain.
%
% inputs:
% -solution: structure containing the following fields
%   -problem: a problem structure containing problem/geometry information
%   -q : density function, vector of size #npts by 2
%
% outputs;

% -[sxx, sxy, syx, syy]: components of stress tensor, evaluated at the 
%quadrature points on the surface

%disp('Evaluating stress...');
local_indices = solution_local.local_indices;

domain = solution.problem.domain;
if solution.problem.periodic
    Lx = domain.Lx;
    Ly = domain.Ly;
end

xsrc = real(solution.problem.domain.z);
ysrc = imag(solution.problem.domain.z);
xtar = xsrc(local_indices);
ytar = ysrc(local_indices);

ntar = length(xtar);

n1 = real(-1i*solution.problem.domain.zp)./abs(solution.problem.domain.zp);
n2 = imag(-1i*solution.problem.domain.zp)./abs(solution.problem.domain.zp);
wazp = solution.problem.domain.wazp;
q = solution.q;

b1 = ones(size(solution_local.problem.domain.z));
b2 = 1i*ones(size(solution_local.problem.domain.z));

if solution.problem.periodic
    % evaluate solution using Ewald, to account for periodic replicates
    [sigmaxx_slp, sigmayx_slp] = StokesSLP_stress_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
        q(:,1).*wazp, q(:,2).*wazp, real(b1), imag(b1), Lx, Ly, 'verbose', 0);

    [sigmaxy_slp, sigmayy_slp] = StokesSLP_stress_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
        q(:,1).*wazp, q(:,2).*wazp, real(b2), imag(b2), Lx, Ly, 'verbose', 0);

    sigma1_slp = sigmaxx_slp + 1i*sigmayx_slp;
    sigma2_slp = sigmaxy_slp + 1i*sigmayy_slp;
    
    % correct with special quadrature
%     sigma1_slp_corrected = stress_slp_on_surface_correction(sigma1_slp, b1, solution_local, type);
%     sigma2_slp_corrected = stress_slp_on_surface_correction(sigma2_slp, b2, solution_local, type);

    [sigma1_slp_corrected,~] = mex_SQ_slp_stress(xtar(:)+1i*(ytar(:)+1e-60),...
        xsrc(:)+1i*(ysrc(:)+1e-60), domain.zp, domain.quad_weights, domain.panel_breaks,...
        domain.wazp, domain.z32, domain.zp32, domain.quad_weights32, domain.wazp32, ...
        solution.q(:,1)+1i*(solution.q(:,2)+1e-60),ones(ntar,1)+1e-60*1i*ones(ntar,1),...
        sigma1_slp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
        domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
        domain.reference_cell,true);
    
    [sigma2_slp_corrected,~] = mex_SQ_slp_stress(xtar(:)+1i*(ytar(:)+1e-60),...
        xsrc(:)+1i*(ysrc(:)+1e-60), domain.zp, domain.quad_weights, domain.panel_breaks,...
        domain.wazp, domain.z32, domain.zp32, domain.quad_weights32, domain.wazp32, ...
        solution.q(:,1)+1i*(solution.q(:,2)+1e-60),1e-60*ones(ntar,1)+1i*ones(ntar,1),...
        sigma2_slp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
        domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
        domain.reference_cell,true);
    
    if ~isinf(solution.problem.eta)
        % add on double-layer potential
        [sigmaxx_dlp, sigmayx_dlp] = StokesDLP_stress_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
            n1, n2, q(:,1).*wazp, q(:,2).*wazp, real(b1), imag(b1), Lx, Ly, 'verbose', 0);

        [sigmaxy_dlp, sigmayy_dlp] = StokesDLP_stress_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
            n1, n2, q(:,1).*wazp, q(:,2).*wazp, real(b2), imag(b2), Lx, Ly, 'verbose', 0);

        sigma1_dlp = sigmaxx_dlp + 1i*sigmayx_dlp;
        sigma2_dlp = sigmaxy_dlp + 1i*sigmayy_dlp;
    
        % correct with special quadrature
%         sigma1_dlp_corrected = stress_dlp_on_surface_correction(sigma1_dlp, b1, solution_local, type);
%         sigma2_dlp_corrected = stress_dlp_on_surface_correction(sigma2_dlp, b2, solution_local, type);
        
        [sigma1_dlp_corrected,~] = mex_SQ_dlp_stress(xtar(:)+1i*(ytar(:)+1e-60),...
            xsrc(:)+1i*(ysrc(:)+1e-60), domain.zp, domain.quad_weights, domain.panel_breaks,...
            domain.wazp, domain.z32, domain.zp32, domain.quad_weights32, domain.wazp32, ...
            solution.q(:,1)+1i*(solution.q(:,2)+1e-60),ones(ntar,1)+1e-60*1i*ones(ntar,1),...
            sigma1_dlp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
            domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
            domain.reference_cell,true);

        [sigma2_dlp_corrected,~] = mex_SQ_dlp_stress(xtar(:)+1i*(ytar(:)+1e-60),...
            xsrc(:)+1i*(ysrc(:)+1e-60), domain.zp, domain.quad_weights, domain.panel_breaks,...
            domain.wazp, domain.z32, domain.zp32, domain.quad_weights32, domain.wazp32, ...
            solution.q(:,1)+1i*(solution.q(:,2)+1e-60),1e-60*ones(ntar,1)+1i*ones(ntar,1),...
            sigma2_dlp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
            domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
            domain.reference_cell,true);
        
        sigma1 = solution.problem.eta*sigma1_slp_corrected + sigma1_dlp_corrected;
        sigma2 = solution.problem.eta*sigma2_slp_corrected + sigma2_dlp_corrected;
    else
        %sigma1 = sigma1_slp;
        %sigma2 = sigma2_slp;
        sigma1 = sigma1_slp_corrected;
        sigma2 = sigma2_slp_corrected;
    end
else
    % TODO
end

sxx = real(sigma1);
sxy = real(sigma2);
syx = imag(sigma1);
syy = imag(sigma2);

% add on driving pressure
if solution.problem.periodic
    sxx = sxx - solution.problem.pressure_gradient_x*xtar;
    sxx = sxx - solution.problem.pressure_gradient_y*ytar;
    syy = syy - solution.problem.pressure_gradient_x*xtar;
    syy = syy - solution.problem.pressure_gradient_y*ytar;
end
