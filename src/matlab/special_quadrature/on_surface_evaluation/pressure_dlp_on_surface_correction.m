function Ph = pressure_dlp_on_surface_correction(P, solution)
%PRESSURE_DLP_ON_SURFACE_CORRECTION corrects the double-layer 
%pressure for on-surface evaluation. Adds on the jump as the target
%approaches the boundary from the exterior of the domain
%
%inputs:
%-P: uncorrected double-layer pressure evaluated at quadrature points on 
%surface
%-solution: solution structure containing geometry information, as well as
%the density 
%
%output:
%-Ph: corrected pressure

domain = solution.problem.domain;

zsrc = domain.z;
zpsrc = domain.zp;
wsrc = domain.quad_weights;

qsrc = solution.q(:,1) + 1i*solution.q(:,2);
Ph = P;
nw = size(domain.wall_indices,1);

wall_start = 1;

for i = 1:nw
    
    indices = domain.wall_indices(i,1):domain.wall_indices(i,2);
    npan = length(indices)/16;
    
    panel_breaks_z = domain.panel_breaks(wall_start:wall_start + npan - 1);
    
    %subract off self contribution
    for j = indices
        
        indices_tmp = indices;
        indices_tmp(indices==j) = [];
        
        Ph(j) = Ph(j) + imag(sum(qsrc(indices_tmp).*wsrc(indices_tmp).*zpsrc(indices_tmp)...
            ./(zsrc(j) - zsrc(indices_tmp)).^2))/pi;
    end
    
    %add on special quadrature
    Ph(indices) = Ph(indices) - imag(hypersingular_on_surface_evaluation(qsrc(indices), ...
        zsrc(indices), zpsrc(indices), wsrc(indices), panel_breaks_z))/pi;
    
    wall_start = wall_start + npan;
end