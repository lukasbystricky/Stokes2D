function [Px, Py] = evaluate_pressure_gradient_on_surface(solution, solution_local, type)
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
local_indices = solution_local.local_indices;

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
    dpslp = StokesSLP_pressure_grad_ewald_2p(xsrc, ysrc, xtar(:), ytar(:),...
        q(:,1).*wazp, q(:,2).*wazp, Lx, Ly,...
        'verbose', 1)';

    % correct with special quadrature
    dpslp = pressure_gradient_slp_on_surface_correction(dpslp(:,1) + 1i*dpslp(:,2), solution_local, type);
    
    if isinf(solution.problem.eta)
        dP = dpslp;
    else
        
        dpdlp = StokesDLP_pressure_grad_ewald_2p(xsrc, ysrc, xtar(:), ytar(:), n1, n2,...
            q(:,1).*wazp, q(:,2).*wazp, Lx, Ly,...
            'verbose', 1)';
        
        % correct with special quadrature
        dpdlp = pressure_gradient_dlp_on_surface_correction(dpdlp(:,1) + 1i*dpdlp(:,2), solution_local, type);
    
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
