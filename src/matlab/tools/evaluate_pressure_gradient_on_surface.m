function [Px, Py] = evaluate_pressure_gradient_on_surface(solution, bodies)
%EVALUTATE_PRESSURE_GRADIENT_ON_SURFACE evaluates the pressure gradient at 
%the quadrature points on the surface of a domain. Adds on the jump 
%corresponding to approaching the boundary from the fluid part of the domain.
%
% inputs:
% -solution: structure containing the following fields
%   -problem: a problem structure containing problem/geometry information
%   -q : density function, vector of size #npts by 2
%
% outputs;
% -dP: pressure gradient, evaluated at the quadrature points on the surface

disp('Evaluating pressure gradient...');

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

n1 = real(-1i*solution.problem.domain.zp(indices))./abs(solution.problem.domain.zp(indices));
n2 = imag(-1i*solution.problem.domain.zp(indices))./abs(solution.problem.domain.zp(indices));
wazp = solution.problem.domain.wazp(indices);
q = solution.q(indices,:);


if solution.problem.periodic
    
    % evaluate solution using Ewald, to account for periodic replicates
    dpslp = StokesSLP_pressure_grad_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
        q(:,1).*wazp, q(:,2).*wazp, Lx, Ly,...
        'verbose', 1)';

    % correct with special quadrature
    dpslp = pressure_gradient_slp_on_surface_correction(dpslp(:,1) + 1i*dpslp(:,2), solution, bodies);
    
    if isinf(solution.problem.eta)
        dP = dpslp;
    else
        
        dpdlp = StokesDLP_pressure_grad_ewald_2p(xsrc, ysrc, xtar(:), ytar(:), n1, n2,...
            q(:,1).*wazp, q(:,2).*wazp, Lx, Ly,...
            'verbose', 1)';
        
        % correct with special quadrature
        dpdlp = pressure_gradient_dlp_on_surface_correction(dpdlp(:,1) + 1i*dpdlp(:,2), solution, bodies);
    
        dP = dpdlp + solution.problem.eta*dpslp;
    end
end

% add on driving pressure
if solution.problem.periodic
    dP = dP + solution.problem.pressure_gradient_x;
    dP = dP + 1i*solution.problem.pressure_gradient_y;
end

Px = real(dP);
Py = imag(dP);
