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

%disp('Evaluating velocity...');

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

if solution.problem.periodic
    % evaluate solution using Ewald, to account for periodic replicates
    [u1_slp, u2_slp] = StokesSLP_ewald_2p(xsrc, ysrc, xtar, ytar, real(qwazp), imag(qwazp),...
        Lx, Ly, 'verbose', 0);

    % correct with special quadrature
    G = u1_slp + 1i*u2_slp;

    wall_indices = domain_local.wall_indices;
    % apply log-quadrature corrections for the Stokeslet
    for i = 1:size(domain_local.wall_indices,1)

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
        [u1_dlp, u2_dlp] = StokesDLP_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, real(qwazp), imag(qwazp),...
            Lx, Ly, 'verbose', 0);

        T = u1_dlp + 1i*u2_dlp;

        % diagonal elements of the double-layer
        u_dlp = T + domain_local.quad_weights.*(q(local_indices).*imag(zpp./zp) + ...
            conj(q(local_indices)).*imag(zpp.*conj(zp))./conj(zp).^2)/(4*pi);

        u = solution.problem.eta*u_slp + u_dlp - q(local_indices)/2;
    else
        u = u_slp;
    end
else % completion formulation
    fmm = 0;
    q1 = solution.q(:,1);
    q2 = solution.q(:,2);
    wazp = domain_local.wazp;
    if fmm
        q1tmp = q1.*wazp;
        q2tmp = q2.*wazp;
        [u1_dlp,u2_dlp] = stokesDLPfmm(q1tmp(:),q2tmp(:),xsrc(:),ysrc(:),n1(:),n2(:));
        
        % note negative sign in front of double-layer
        u_dlp = -u1_dlp - 1i*u2_dlp;
    else
        u1_dlp = zeros(numel(xtar),1);
        u2_dlp = zeros(numel(ytar),1);
        indices = 1:numel(xtar);
        for k = 1:numel(xtar)
            
            % skip self-interaction
            indices_tmp = indices;
            indices_tmp(indices==k) = [];
            
            rx = xtar(k) - xsrc(indices_tmp);
            ry = ytar(k) - ysrc(indices_tmp);
            rho4 = (rx.^2 + ry.^2).^2;
            
            rdotq = rx.*q1(indices_tmp).*wazp(indices_tmp) + ...
                ry.*q2(indices_tmp).*wazp(indices_tmp);
            rdotn = rx.*n1(indices_tmp) + ry.*n2(indices_tmp);
            
            u1_dlp(k) = 4*sum(rdotn.*rdotq./rho4.*rx);
            u2_dlp(k) = 4*sum(rdotn.*rdotq./rho4.*ry);

        end
        
        u1_dlp = u1_dlp/4/pi;
        u2_dlp = u2_dlp/4/pi;
        
        u_dlp = -u1_dlp - 1i*u2_dlp;
        
    end

    % diagonal elements of the double-layer
    % Not sure why we have to split into the two cases. DK
    u = u_dlp;
    outer_wall_indices = domain.wall_indices(1,1):domain.wall_indices(1,2);
    inner_wall_indices = domain.wall_indices(2,1):domain.wall_indices(end,2);
    quad_weights_in = domain.quad_weights(inner_wall_indices);
    zp_in = domain.zp(inner_wall_indices);
    zpp_in = domain.zpp(inner_wall_indices);
    q_in = q(inner_wall_indices);
    quad_weights_out = domain.quad_weights(outer_wall_indices);
    zp_out = domain.zp(outer_wall_indices);
    zpp_out = domain.zpp(outer_wall_indices);
    q_out = q(outer_wall_indices);
    u(inner_wall_indices) = u(inner_wall_indices) + quad_weights_in.*...
        (q_in.*imag(zpp_in./zp_in) + ...
        conj(q_in).*imag(zpp_in.*conj(zp_in))./conj(zp_in).^2)/(4*pi);
    u(outer_wall_indices) = u(outer_wall_indices) - quad_weights_out.*...
        (q_out.*imag(zpp_out./zp_out) + ...
        conj(q_out).*imag(zpp_out.*conj(zp_out))./conj(zp_out).^2)/(4*pi);

    % remove container nullspace
    % Not sure if needed. DK
%     outer_wall_indices = domain.wall_indices(1,1):domain.wall_indices(1,2);
%     n_outer1 = n1(outer_wall_indices);
%     n_outer2 = n2(outer_wall_indices);
%     qw_outer1 = real(qwazp(outer_wall_indices));
%     qw_outer2 = imag(qwazp(outer_wall_indices));
%     N0 = (n_outer1 + 1i*n_outer2).*(n_outer1.*qw_outer1 + n_outer2.*qw_outer2);
%     u(outer_wall_indices) = u(outer_wall_indices) + N0;
    
    % add on completion flow from rotlets and Stokeslets
    [uS, uR] = completion_contribution(domain.centers(2:end), xtar+1i*ytar,...
        solution.forces, solution.torques);

    u = -q(local_indices)/2 + u + uS + uR;

end

U = real(u);
V = imag(u);

if solution.problem.periodic
    U = U + solution.u_avg(1);
    V = V + solution.u_avg(2);
end
