function [U, V] = evaluate_velocity_on_surface(solution, bodies)
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
if solution.problem.periodic
    Lx = domain.Lx;
    Ly = domain.Ly;
end

indices = [];
for i = bodies
    indices = [indices, solution.problem.domain.wall_indices(bodies,1) : ...
                solution.problem.domain.wall_indices(bodies,2)];
end

xsrc = real(solution.problem.domain.z(indices));
ysrc = imag(solution.problem.domain.z(indices));
xtar = xsrc;
ytar = ysrc;
zp = domain.zp(indices);
zpp = domain.zpp(indices);

n1 = real(-1i*solution.problem.domain.zp(indices))./abs(solution.problem.domain.zp(indices));
n2 = imag(-1i*solution.problem.domain.zp(indices))./abs(solution.problem.domain.zp(indices));
wazp = solution.problem.domain.wazp(indices);
q = solution.q(indices,1)+1i*solution.q(indices,2);
qwazp = q.*wazp;

% evaluate solution using Ewald, to account for periodic replicates
[u_slp1, u_slp2] = StokesSLP_ewald_2p(xsrc, ysrc, xtar, ytar, real(qwazp), imag(qwazp),...
    Lx, Ly);

% correct with special quadrature
G = u_slp1 + 1i*u_slp2;

% apply log-quadrature corrections for the Stokeslet
wall_indices = domain.wall_indices;
for i = 1:bodies
    
    panels_per_wall = (wall_indices(i,2)-wall_indices(i,1)+1)/16;
    thiswall = wall_indices(i,1):wall_indices(i,2);
    G(thiswall) = G(thiswall) - qwazp(thiswall).*...
        log(pi/panels_per_wall*abs(domain.zp(thiswall)))/(4*pi);
    
    for j = 1:panels_per_wall
        
        source_ind = wall_indices(i,1)-1+((j-1)*16+(1:16));
        target_ind = wall_indices(i,1)+...
            mod((j-1)*16+(-3:20)+panels_per_wall*16-1, panels_per_wall*16);
        
        qwazp_tmp = real(solution.q(source_ind,1) + 1i*solution.q(source_ind,2)).*domain.wazp(source_ind);
        correction = -domain.Lmod*(qwazp_tmp) / (4*pi);
        
        G(target_ind) = G(target_ind) + correction;
    end
end

u_slp = G + (q + zp./conj(zp).*conj(q)).*wazp/(8*pi);

if ~isinf(solution.problem.eta)
    % add on double-layer potential
    [u_dlp1, u_dlp2] = StokesDLP_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, real(qwazp), imag(qwazp),...
        Lx, Ly);
    
    T = u_dlp1 + 1i*u_dlp2;
    
    % diagonal elements of the double-layer
    u_dlp = T + domain.quad_weights(indices).*(q.*imag(zpp./zp) + ...
        conj(q).*imag(zpp.*conj(zp))./conj(zp).^2)/(4*pi);
    
    u = solution.problem.eta*u_slp + u_dlp - q/2;
else
    u = u_slp;
end

U = real(u) + solution.u_avg(1);
V = imag(u) + solution.u_avg(2);