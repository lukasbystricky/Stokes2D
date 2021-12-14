function omegac = vorticity_slp_on_surface_correction(omega, solution_local, type)
%VORTICITY_GRADIENT_SLP_ON_SURFACE_CORRECTION corrects the single-layer 
%vorticity for on-surface evaluation. Adds on the jump as the 
%target approaches the boundary from the fluid part of the domain.
%
%inputs:
%-omega: single-layer vorticity evaluated at quadrature points on 
%surface
%-solution: solution structure containing geometry information, as well as
%the density 
%
%output:
%-omegac: corrected vorticity

domain = solution_local.problem.domain;

zsrc = domain.z;
zpsrc = domain.zp;
wsrc = domain.quad_weights;

nsrc = -1i*zpsrc./abs(zpsrc);
qsrc = solution_local.q(:,1) + 1i*solution_local.q(:,2);

omegac = omega;

wall_start = 1;
for i = 1:size(domain.wall_indices,1)
    
    indices = domain.wall_indices(i,1):domain.wall_indices(i,2);
    npan = length(indices)/16;
    
    panel_breaks_z = domain.panel_breaks(wall_start:wall_start + npan - 1);
    
    %subract off self contribution
    for j = indices
        
        indices_tmp = indices;
        indices_tmp(indices==j) = [];
        
        omegac(j) = omegac(j) + real(sum(qsrc(indices_tmp).*wsrc(indices_tmp).*...
            zpsrc(indices_tmp)./(nsrc(indices_tmp).*(zsrc(indices_tmp) - zsrc(j)))))/(2*pi);
    end
    
    %add on special quadrature
    % different sign due to implementation of the Cauchy integral
    omegac(indices) = omegac(indices) + real(cauchy_on_surface_evaluation(qsrc(indices)./nsrc(indices),...
        zsrc(indices),zpsrc(indices),wsrc(indices),panel_breaks_z,type))/(2*pi);
     
   wall_start = wall_start + npan + 1;       
end