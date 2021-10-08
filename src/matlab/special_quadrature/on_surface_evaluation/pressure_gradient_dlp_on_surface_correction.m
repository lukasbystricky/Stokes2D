function dPc = pressure_gradient_dlp_on_surface_correction(dP, solution_local, type)
%PRESSURE_GRADIENT_DLP_ON_SURFACE_CORRECTION corrects the double-layer 
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
%-dPc: corrected pressure gradient

domain = solution_local.problem.domain;

zsrc = domain.z;
zpsrc = domain.zp;
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
        
        dPc(j) = dPc(j) - 2*1i*conj(sum(qsrc(indices_tmp).*wsrc(indices_tmp).*...
                zpsrc(indices_tmp)./(zsrc(j) - zsrc(indices_tmp)).^3))/pi;
    end    
    
    %add on special quadrature
    dPc(indices) = dPc(indices) + 2*1i*conj(supersingular_on_surface_evaluation(qsrc(indices), ...
                zsrc(indices), zpsrc(indices), wsrc(indices), panel_breaks_z, type))/pi; 
            
   wall_start = wall_start + npan + 1;       
end