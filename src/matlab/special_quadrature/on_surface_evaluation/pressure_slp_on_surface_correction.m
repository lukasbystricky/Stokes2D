function Pc = pressure_slp_on_surface_correction(P, solution_local, type)
%PRESSURE_SLP_ON_SURFACE_CORRECTION corrects the single-layer 
%pressure for on-surface evaluation. Adds on the jump as the target
%approaches the boundary from the fluid part of the domain.
%
%inputs:
%-P: uncorrected single-layer pressure evaluated at quadrature points on 
%surface
%-solution: solution structure containing geometry information, as well as
%the density 
%
%output:
%-Pc: corrected pressure

domain = solution_local.problem.domain;

zsrc = domain.z;
zpsrc = domain.zp;
nsrc = -1i*zpsrc./abs(zpsrc);
wsrc = domain.quad_weights;

qsrc = (solution_local.q(:,1) + 1i*solution_local.q(:,2))./nsrc;

Pc = P;

wall_start = 1;
for i = size(domain.wall_indices,1)

    
    indices = domain.wall_indices(i,1):domain.wall_indices(i,2);
    npan = length(indices)/16;
    
    panel_breaks_z = domain.panel_breaks(wall_start:wall_start + npan - 1);
    
    %subract off self contribution
    for j = indices
        
        indices_tmp = indices;
        indices_tmp(indices==j) = [];
        
        Pc(j) = Pc(j) - imag(sum(qsrc(indices_tmp).*wsrc(indices_tmp).*...
                zpsrc(indices_tmp)./(zsrc(j) - zsrc(indices_tmp))))/(2*pi);
    end    
    
    %add on special quadrature
    Pc(indices) = Pc(indices) + imag(cauchy_on_surface_evaluation(qsrc(indices), ...
                zsrc(indices), zpsrc(indices), wsrc(indices), panel_breaks_z, type))/(2*pi); 
            
   wall_start = wall_start + npan + 1;       
end