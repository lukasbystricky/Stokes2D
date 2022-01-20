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

disp('Evaluating stress...');
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

n1 = real(-1i*solution.problem.domain.zp)./abs(solution.problem.domain.zp);
n2 = imag(-1i*solution.problem.domain.zp)./abs(solution.problem.domain.zp);
wazp = solution.problem.domain.wazp;
q = solution.q;

b1 = ones(size(solution_local.problem.domain.z));
b2 = 1i*ones(size(solution_local.problem.domain.z));

if solution.problem.periodic
    % evaluate solution using Ewald, to account for periodic replicates
    [sigmaxx_slp, sigmayx_slp] = StokesSLP_stress_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
        q(:,1).*wazp, q(:,2).*wazp, real(b1), imag(b1), Lx, Ly, 'verbose', 1);

    [sigmaxy_slp, sigmayy_slp] = StokesSLP_stress_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
        q(:,1).*wazp, q(:,2).*wazp, real(b2), imag(b2), Lx, Ly, 'verbose', 1);

    sigma_slp1 = sigmaxx_slp + 1i*sigmayx_slp;
    sigma_slp2 = sigmaxy_slp + 1i*sigmayy_slp;
    
    % correct with special quadrature
    sigma_slp1 = stress_slp_on_surface_correction(sigma_slp1, b1, solution_local, type);
    sigma_slp2 = stress_slp_on_surface_correction(sigma_slp2, b2, solution_local, type);

    if ~isinf(solution.problem.eta)
        % add on double-layer potential
        [sigmaxx_dlp, sigmayx_dlp] = StokesDLP_stress_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
            n1, n2, q(:,1).*wazp, q(:,2).*wazp, real(b1), imag(b1), Lx, Ly, 'verbose', 1);

        [sigmaxy_dlp, sigmayy_dlp] = StokesDLP_stress_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
            n1, n2, q(:,1).*wazp, q(:,2).*wazp, real(b2), imag(b2), Lx, Ly, 'verbose', 1);

        sigma_dlp1 = sigmaxx_dlp + 1i*sigmayx_dlp;
        sigma_dlp2 = sigmaxy_dlp + 1i*sigmayy_dlp;
    
        % correct with special quadrature
        sigma_dlp1 = stress_dlp_on_surface_correction(sigma_dlp1, b1, solution_local, type);
        sigma_dlp2 = stress_dlp_on_surface_correction(sigma_dlp2, b2, solution_local, type);

        sigma1 = solution.problem.eta*sigma_slp1 + sigma_dlp1;
        sigma2 = solution.problem.eta*sigma_slp2 + sigma_dlp2;
    else
        sigma1 = sigma_slp1;
        sigma2 = sigma_slp2;
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
