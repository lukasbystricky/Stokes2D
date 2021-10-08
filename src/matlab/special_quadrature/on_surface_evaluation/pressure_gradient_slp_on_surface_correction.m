function dPc = pressure_gradient_slp_on_surface_correction(dP, solution_local, type)
%PRESSURE_GRADIENT_SLP_ON_SURFACE_CORRECTION corrects the single-layer 
%pressure gradient for on-surface evaluation. Adds on the jump as the 
%target approaches the boundary from the fluid part of the domain.
%
%inputs:
%-dP: uncorrected single-layer pressure evaluated at quadrature points on 
%surface
%-solution: solution structure containing geometry information, as well as
%the density 
%
%output:
%-dPc: corrected pressure

domain = solution_local.problem.domain;

zsrc = domain.z;
zpsrc = domain.zp;
nsrc = -1i*zpsrc./abs(zpsrc);
wsrc = domain.quad_weights;

qsrc = (solution_local.q(:,1) + 1i*solution_local.q(:,2));

dPc = dP;

wall_start = 1;
for i = 1:size(domain.wall_indices,1)
    
    indices = domain.wall_indices(i,1):domain.wall_indices(i,2);
    npan = length(indices)/16;
    
    panel_breaks_z = domain.panel_breaks(wall_start:wall_start + npan - 1);
    
    %subract off self contribution
    for j = indices
        
        indices_tmp = indices;
        indices_tmp(indices==j) = [];
        
        dPc(j) = dPc(j) + 1i*conj(sum(qsrc(indices_tmp).*wsrc(indices_tmp).*...
                zpsrc(indices_tmp)./(nsrc(indices_tmp).*(zsrc(j) - zsrc(indices_tmp)).^2)))/(2*pi);
    end    
    
    %add on special quadrature
    dPc(indices) = dPc(indices) - 1i*conj(hypersingular_on_surface_evaluation(qsrc(indices)./nsrc(indices), ...
                zsrc(indices), zpsrc(indices), wsrc(indices), panel_breaks_z, type))/(2*pi); 
            
   wall_start = wall_start + npan + 1;       
end