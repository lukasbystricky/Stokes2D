function [P, xtar, ytar] = evaluate_pressure_on_surface(solution, solution_local, type)
%EVALUTATE_PRESSURE_ON_SURFACE evaluates the pressure at the quadrature
%points on the surface of a domain. Adds on the jump corresponding to
%approaching the boundary from the fluid part of the domain.
%
% inputs:
% -solution: structure containing the following fields
%   -problem: a problem structure containing problem/geometry information
%   -q : density function, vector of size #npts by 2
%
% outputs;
% -P: pressure, evaluated at the quadrature points on the surface

disp('Evaluating pressure...');
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


if solution.problem.periodic
    
    % evaluate solution using Ewald, to account for periodic replicates
    pslp = StokesSLP_pressure_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
        q(:,1).*wazp, q(:,2).*wazp, Lx, Ly,...
        'verbose', 1)';

    pslp = pressure_slp_on_surface_correction(pslp, solution_local, type);

    if isinf(solution.problem.eta)
        P = pslp;
    else
        
        pdlp = StokesDLP_pressure_ewald_2p(xsrc, ysrc, xtar(:), ytar(:), n1, n2,...
            q(:,1).*wazp, q(:,2).*wazp, Lx, Ly,...
            'verbose', 1)';
        
        % correct with special quadrature
        pdlp = pressure_dlp_on_surface_correction(pdlp, solution_local, type);
    
        P = pdlp + solution.problem.eta*pslp;
    end
    
else % not periodic, only double-layer plus Stokeslets
    
    % note that the pressure is evaluated without an FMM, so it will be
    % quite slow to evaluate if the number of target points is large
    pdlp = evaluate_double_layer_pressure_direct(xtar(:), ytar(:), domain.z, domain.zp,...
        solution.q(:,1)+1i*solution.q(:,2), domain.quad_weights);
    
    % correct with special quadrature
    % TODO
    
    pstokeslets = zeros(size(xtar(:)));
    for i = 1:length(solution.forces)
        R = (xtar(:) - domain.centers(i,1)) + 1i*(ytar(:) - domain.centers(i,2));
        pstokeslets = pstokeslets + real(conj(solution.forces(i))*R)./abs(R).^2;
    end
    
    P = pdlp + pstokeslets;
    
end

% add on driving pressure
if solution.problem.periodic
    P = P + solution.problem.pressure_gradient_x*xtar;
    P = P + solution.problem.pressure_gradient_y*ytar;
end