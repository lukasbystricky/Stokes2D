function [Ux, Uy, Vx, Vy] = evaluate_velocity_gradient_on_surface(solution, solution_local, type)
%EVALUTATE_VELOCITY_GRADIENT_ON_SURFACE evaluates the velocity gradient at 
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
% -[Ux, Uy, Vx, Vy]: components of the velocity gradient, evaluated at the 
%quadrature points on the surface

disp('Evaluating velocity gradient...');
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

% evaluate solution using Ewald, to account for periodic replicates
b1 = ones(size(solution_local.problem.domain.z));
b2 = 1i*ones(size(solution_local.problem.domain.z));

[ux, vx] = StokesSLP_gradient_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
   q(:,1).*wazp, q(:,2).*wazp, real(b1), imag(b1), ...
    Lx, Ly,'verbose', 1);

[uy, vy] = StokesSLP_gradient_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
    q(:,1).*wazp, q(:,2).*wazp, real(b2), imag(b2), ...
    Lx, Ly,'verbose', 1);

% correct with special quadrature
ugrad_slp1 = velocity_gradient_slp_on_surface_correction(ux +...
                1i*vx, b1, solution_local, type);
ugrad_slp2 = velocity_gradient_slp_on_surface_correction(uy +...
                1i*vy, b2, solution_local, type);

if ~isinf(solution.problem.eta)
    % add on double-layer potential
    
    [ux, vx] = StokesDLP_gradient_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
        n1, n2, q(:,1).*wazp, q(:,2).*wazp, real(b1), imag(b1), ...
        Lx, Ly,'verbose', 1);
    
    [uy, vy] = StokesDLP_gradient_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
        n1, n2, q(:,1).*wazp, q(:,2).*wazp, real(b2), imag(b2), ...
        Lx, Ly,'verbose', 1);
    
    % correct with special quadrature
    ugrad_dlp1 = velocity_gradient_dlp_on_surface_correction(ux + 1i*vx, b1, solution_local, type);
    ugrad_dlp2 = velocity_gradient_dlp_on_surface_correction(uy + 1i*vy, b2, solution_local, type);
    
    ugrad1 = solution.problem.eta*ugrad_slp1 + ugrad_dlp1;
    ugrad2 = solution.problem.eta*ugrad_slp2 + ugrad_dlp2;
else
    ugrad1 = ugrad_slp1;
    ugrad2 = ugrad_slp2;
end


Ux = real(ugrad1);
Uy = real(ugrad2);
Vx = imag(ugrad1);
Vy = imag(ugrad2);
