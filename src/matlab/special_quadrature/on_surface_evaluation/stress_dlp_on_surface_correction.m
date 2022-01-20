function sigmac = stress_dlp_on_surface_correction(sigma, btar, solution_local, type)
%STRESS_DLP_ON_SURFACE_CORRECTION corrects the double-layer 
%stress for on-surface evaluation. Returns a vector which is the
%gradient dotted with a direction vector b. Adds on the jump as the target
%approaches the boundary from the fluid part of the domain.
%
%inputs:
%-sigma: uncorrected double-layer stress evaluated at quadrature 
%points on surface
%-btar: direction vector at every target point (e.g. normal vector)
%-solution: solution structure containing geometry information, as well as
%the density 
%
%output:
%-sigmac: corrected stress

domain = solution_local.problem.domain;

zsrc = domain.z;
zpsrc = domain.zp;

wsrc = domain.quad_weights;

nsrc = -1i*zpsrc./abs(zpsrc);
qsrc = solution_local.q(:,1) + 1i*solution_local.q(:,2);
sigmac = sigma;

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
        
        Ih1 = sum((qtmp.*wtmp.*zptmp)./r.^2);
        Ih2 = sum((real(qtmp.*conj(ntmp))./(ntmp.*r.^2)).*wtmp.*zptmp);
        Is1 = -sum((conj(ztmp).*qtmp./r.^3).*wtmp.*zptmp);
        Is2 = -sum((qtmp./r.^3).*wtmp.*zptmp);
        
        sigmac(j) = sigmac(j) - (-(btar(j)*imag(Ih1)/2 - 1i*conj(btar(j)*(Ih2+Is1-conj(zsrc(j))*Is2)))/pi);

    end
    
    %add on special quadrature
    ztmp = zsrc(indices);
    qtmp = qsrc(indices);
    ntmp = nsrc(indices);
    zptmp = zpsrc(indices);
    wtmp = wsrc(indices);
    btmp = btar(indices);
    
    % minus sign on super singular due to how it is implemented
    Ih1 = hypersingular_on_surface_evaluation(qtmp,ztmp,zptmp,wtmp,panel_breaks_z,type);
    Ih2 = hypersingular_on_surface_evaluation(real(qtmp.*conj(ntmp))./ntmp,ztmp,zptmp,wtmp,panel_breaks_z,type);
    Is1 = -supersingular_on_surface_evaluation(conj(ztmp).*qtmp,ztmp,zptmp,wtmp,panel_breaks_z,type);
    Is2 = -supersingular_on_surface_evaluation(qtmp,ztmp,zptmp,wtmp,panel_breaks_z,type);
    
    sigmac(indices) = sigmac(indices) + (-(btmp.*imag(Ih1)/2 - 1i*conj(btmp.*(Ih2+Is1-conj(ztmp).*Is2)))/pi);
     
    wall_start = wall_start + npan + 1;
end
