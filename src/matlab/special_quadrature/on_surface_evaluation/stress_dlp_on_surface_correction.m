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
periodic = solution_local.problem.periodic;

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
    
    panel_breaks_z = domain.panel_breaks(wall_start:wall_start + npan);

    % check if boundary is a closed curve
    if abs(panel_breaks_z(1)-panel_breaks_z(end)) < 1e-12
        closed_curve = 1;
    else
        closed_curve = 0;
    end
    
    % if problem is periodic and the domain is not a closed curve, then we
    % need to handle periodic replicates in a specific way
    periodic_rep = periodic && ~closed_curve;

    % subract off contribution from the current wall, compensates for 
    % periodic replicates due to them being included in Ewald, which will
    % be computed later using a special quadrature rule
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

        if periodic_rep
            % we also correct for the adjacent panels that have been
            % periodically replicated later, so these need to be removed too
            indices_tmp = [];

            if (j - (i-1)*16*npan) <= 16
                indices_tmp = indices(end-15:end);

                if real(zsrc(indices(2))) > real(zsrc(indices(1)))
                    ztmp = zsrc(indices_tmp) - solution_local.problem.Lx;
                else
                    ztmp = zsrc(indices_tmp) + solution_local.problem.Lx;
                end

            elseif (j - (i-1)*16*npan) >= (length(indices) - 15)
               indices_tmp = indices(1:16);

                if real(zsrc(indices(2))) > real(zsrc(indices(1)))
                    ztmp = zsrc(indices_tmp) + solution_local.problem.Lx;
                else
                    ztmp = zsrc(indices_tmp) - solution_local.problem.Lx;
               end
            end
            
            if ~isempty(indices_tmp)
                r = zsrc(j) - ztmp;
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
        end
    end
    
    %add on special quadrature
    ztmp = zsrc(indices);
    qtmp = qsrc(indices);
    ntmp = nsrc(indices);
    zptmp = zpsrc(indices);
    wtmp = wsrc(indices);
    btmp = btar(indices);
    
    % minus sign on super singular due to how it is implemented
    Ih1 = hypersingular_on_surface_evaluation(qtmp,ztmp,zptmp,wtmp,panel_breaks_z,type,periodic_rep);
    Ih2 = hypersingular_on_surface_evaluation(real(qtmp.*conj(ntmp))./ntmp,ztmp,zptmp,wtmp,panel_breaks_z,type,periodic_rep);
    Is1 = -supersingular_on_surface_evaluation_zsrc(qtmp,ztmp,zptmp,wtmp,panel_breaks_z,type,periodic_rep);
    Is2 = -supersingular_on_surface_evaluation(qtmp,ztmp,zptmp,wtmp,panel_breaks_z,type,periodic_rep);
    
    sigmac(indices) = sigmac(indices) + (-(btmp.*imag(Ih1)/2 - 1i*conj(btmp.*(Ih2+Is1-conj(ztmp).*Is2)))/pi);
     
    wall_start = wall_start + npan + 1;
end
