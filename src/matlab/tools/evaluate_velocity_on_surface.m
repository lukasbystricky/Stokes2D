function [U, V] = evaluate_velocity_on_surface(solution, solution_local)
%EVALUTATE_VELOCITY_ON_SURFACE evaluates the velocity at the quadrature
%points on the surface of a domain. Adds on the jump corresponding to
%approaching the boundary from the fluid part of the domain.
%
% inputs:
% -solution: structure containing the following fields
%   -problem: a problem structure containing problem/geometry information
%   -q : density function, vector of size #npts by 2
%
% outputs;
% -[U,V]: velocity components, evaluated at the quadrature points on the 
%surface

disp('Evaluating velocity...');

domain = solution.problem.domain;
domain_local = solution_local.problem.domain;
if solution.problem.periodic
    Lx = domain.Lx;
    Ly = domain.Ly;
end

local_indices = solution_local.local_indices;

xsrc = real(solution.problem.domain.z);
ysrc = imag(solution.problem.domain.z);
xtar = xsrc(local_indices);
ytar = ysrc(local_indices);
zp = domain.zp(local_indices);
zpp = domain.zpp(local_indices);

n1 = real(-1i*solution.problem.domain.zp)./abs(solution.problem.domain.zp);
n2 = imag(-1i*solution.problem.domain.zp)./abs(solution.problem.domain.zp);
wazp = solution.problem.domain.wazp;
q = solution.q(:,1)+1i*solution.q(:,2);
qwazp = q.*wazp;

qwazp_local = qwazp(local_indices);

% evaluate solution using Ewald, to account for periodic replicates
[u_slp1, u_slp2] = StokesSLP_ewald_2p(xsrc, ysrc, xtar, ytar, real(qwazp), imag(qwazp),...
    Lx, Ly);

% correct with special quadrature
G = u_slp1 + 1i*u_slp2;

wall_indices = domain_local.wall_indices;
% apply log-quadrature corrections for the Stokeslet
for i = size(domain_local.wall_indices,1)
    
    thiswall = wall_indices(i,1):wall_indices(i,2);
    panels_per_wall = (wall_indices(i,2)-wall_indices(i,1)+1)/16;
    
    G(thiswall) = G(thiswall) - qwazp_local(thiswall).*...
        log(pi/panels_per_wall*abs(domain_local.zp(thiswall)))/(4*pi);
    
    for j = 1:panels_per_wall
        
        source_ind = wall_indices(i,1)-1+((j-1)*16+(1:16));
        target_ind = wall_indices(i,1)+...
            mod((j-1)*16+(-3:20)+panels_per_wall*16-1, panels_per_wall*16);
        
        qwazp_tmp = (solution_local.q(source_ind,1) + ...
                    1i*solution_local.q(source_ind,2)).*domain_local.wazp(source_ind);
        correction = -domain.Lmod*qwazp_tmp / (4*pi);
        
        G(target_ind) = G(target_ind) + correction;
    end
end

u_slp = G + (q(local_indices) + zp./conj(zp).*conj(q(local_indices))).*wazp(local_indices)/(8*pi);

if ~isinf(solution.problem.eta)
    % add on double-layer potential
    [u_dlp1, u_dlp2] = StokesDLP_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, real(qwazp), imag(qwazp),...
        Lx, Ly);
    
    T = u_dlp1 + 1i*u_dlp2;
    
    % diagonal elements of the double-layer
    u_dlp = T + domain_local.quad_weights.*(q(local_indices).*imag(zpp./zp) + ...
        conj(q(local_indices)).*imag(zpp.*conj(zp))./conj(zp).^2)/(4*pi);
    
    u = solution.problem.eta*u_slp + u_dlp - q(local_indices)/2;
else
    u = u_slp;
end

U = real(u) + solution.u_avg(1);
V = imag(u) + solution.u_avg(2);