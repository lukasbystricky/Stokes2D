function dVh = velocity_gradient_slp_on_surface_correction(dV, btar, solution_local, type)
%VELOCITY_GRADIENT_SLP_ON_SURFACE_CORRECTION corrects the single-layer 
%velocity gradient for on-surface evaluation. Returns a vector which is the
%gradient dotted with a direction vector b. Adds on the jump as the target
%approaches the boundary from the fluid part of the domain.
%
%inputs:
%-dV: uncorrected single-layer velocity gradient evaluated at quadrature 
%points on surface
%-btar: direction vector at every target point (e.g. normal vector)
%-solution: solution structure containing geometry information, as well as
%the density 
%
%output:
%-dVh: corrected velocity gradient

domain = solution_local.problem.domain;

zsrc = domain.z;
zpsrc = domain.zp;
wsrc = domain.quad_weights;

nsrc = -1i*zpsrc./abs(zpsrc);
qsrc = solution_local.q(:,1) + 1i*solution_local.q(:,2);
dVh = dV;

wall_start = 1;

for i = 1:size(domain.wall_indices,1)
    
    indices = domain.wall_indices(i,1):domain.wall_indices(i,2);
    npan = (domain.wall_indices(i,2)-domain.wall_indices(i,1)+1)/16;
    
    panel_breaks_z = domain.panel_breaks(wall_start:wall_start + npan - 1);
    
    %subract off self contribution
    for j = indices
        
        indices_tmp = indices;
        indices_tmp(indices==j) = [];
        
        r = zsrc(j) - zsrc(indices_tmp);
        ztmp = zsrc(indices_tmp);
        qtmp = qsrc(indices_tmp);
        wtmp = wsrc(indices_tmp);
        ntmp = nsrc(indices_tmp);
        zptmp = zpsrc(indices_tmp);
        
        dVh(j) = dVh(j) - (-1i*conj(btar(j))/2*conj(sum((conj(zsrc(j))*qtmp./(ntmp.*r.^2) - ...
            conj(ztmp).*qtmp./(ntmp.*r.^2) + conj(qtmp)./(ntmp.*r)).*wtmp.*zptmp)) + ...
            1i*btar(j)*real(sum(qtmp./(ntmp.*r).*wtmp.*zptmp)))/(4*pi);
    end
    
    %add on special quadrature
    ztmp = zsrc(indices);
    qtmp = qsrc(indices);
    ntmp = nsrc(indices);
    zptmp = zpsrc(indices);
    wtmp = wsrc(indices);
    btmp = btar(indices);

    dVh(indices) = dVh(indices) + (-1i*conj(btmp)/2.*conj(...
         conj(ztmp).*hypersingular_on_surface_evaluation(qtmp./ntmp, ztmp, zptmp, wtmp, panel_breaks_z, type) - ...
         hypersingular_on_surface_evaluation(conj(ztmp).*qtmp./ntmp, ztmp, zptmp, wtmp, panel_breaks_z, type) + ...
         cauchy_on_surface_evaluation(conj(qtmp)./ntmp, ztmp, zptmp, wtmp, panel_breaks_z, type)) + ...
         1i*btmp.*real(cauchy_on_surface_evaluation(qtmp./ntmp, ztmp, zptmp, wtmp, panel_breaks_z, type)))/(4*pi);
    
    wall_start = wall_start + npan + 1;
end