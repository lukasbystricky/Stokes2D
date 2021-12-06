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
omega=0;
% certain things are easier in complex notation
zsrc = xsrc + 1i*ysrc;
ztar = xtar + 1i*ytar;
qc = solution.q(:,1)+1i*solution.q(:,2);
N = length(ztar);

if solution.problem.periodic
    % TODO
else
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