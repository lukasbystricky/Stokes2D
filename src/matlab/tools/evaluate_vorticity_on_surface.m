function [omega, xtar, ytar] = evaluate_vorticity_on_surface(solution, solution_local, type)
%EVALUTATE_VORTICITY_ON_SURFACE evaluates the vorticity at the quadrature
%points on the surface of a domain. Adds on the jump corresponding to
%approaching the boundary from the fluid part of the domain.
%
% inputs:
% -solution: structure containing the following fields
%   -problem: a problem structure containing problem/geometry information
%   -q : density function, vector of size #npts by 2
%
% outputs;
% -omega: vorticity, evaluated at the quadrature points on the surface

disp('Evaluating vorticity...');
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

% certain things are easier in complex notation
zsrc = xsrc + 1i*ysrc;
ztar = xtar + 1i*ytar;
n1 = real(-1i*solution.problem.domain.zp)./abs(solution.problem.domain.zp);
n2 = imag(-1i*solution.problem.domain.zp)./abs(solution.problem.domain.zp);
wazp = solution.problem.domain.wazp;
qc = solution.q(:,1)+1i*solution.q(:,2);
N = length(ztar);

if solution.problem.periodic
    % evaluate solution using Ewald, to account for periodic replicates
    omega_slp = StokesSLP_vorticity_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
        solution.q(:,1).*wazp, solution.q(:,2).*wazp, Lx, Ly,...
        'verbose', 1);
    
    % correct with special quadrature
    %omega_slp_corrected = vorticity_slp_on_surface_correction(omega_slp', solution_local, type);
    omega_slp_corrected = vorticity_slp_on_surface_correction_new(omega_slp', solution_local, type);

    if isinf(solution.problem.eta)
        omega = omega_slp_corrected;
    else
        % add on double-layer potential
        omega_dlp = StokesDLP_vorticity_ewald_2p(xsrc, ysrc, xtar(:), ytar(:), n1, n2,...
            solution.q(:,1).*wazp, solution.q(:,2).*wazp, Lx, Ly,...
            'verbose', 1);
        
        % correct with special quadrature
        %omega_dlp_corrected = vorticity_dlp_on_surface_correction(omega_dlp', solution_local, type);
        omega_dlp_corrected = vorticity_dlp_on_surface_correction_new(omega_dlp', solution_local, type);
        
        omega = omega_dlp_corrected + solution.problem.eta*omega_slp_corrected;
    end
else
    % this has to be looked over
    omega_dlp = zeros(N,1);
    for k = 1:N
        % skip self interaction term
        ind = [(1:k-1) (k+1:N)];
        
        rho = zsrc(ind) - ztar(k);
        omega_dlp(k) = sum(qc(ind).*domain.zp(ind).*domain.quad_weights(ind)./rho.^2);
    end
    omega_dlp = -real(omega_dlp)/pi;
    
    % correct with special quadrature
    omega_dlp_corrected = vorticity_dlp_on_surface_correction(omega_dlp, solution_local, type);

    % add on Stokeslet and rotlet contribution
    omegaR = zeros(N,1);
    omegaS = completion_contribution_vorticity(domain.centers(2:end), ztar,...
        solution.forces);
    
    omega = omega_dlp_corrected + omegaS + omegaR;
end