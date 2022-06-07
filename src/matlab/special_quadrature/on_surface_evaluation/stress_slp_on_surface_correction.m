function sigmac = stress_slp_on_surface_correction(sigma, btar, solution_local, type)
%STRESS_SLP_ON_SURFACE_CORRECTION corrects the single-layer 
%stress for on-surface evaluation. Returns a vector which is the
%gradient dotted with a direction vector b. Adds on the jump as the target
%approaches the boundary from the fluid part of the domain.
%
%inputs:
%-sigma: uncorrected single-layer stress evaluated at quadrature 
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
        
        Ic1 = -sum((qtmp./(ntmp.*r)).*wtmp.*zptmp);
        Ic2 = -sum((conj(qtmp)./(ntmp.*r)).*wtmp.*zptmp);
        Ih1 = sum((conj(ztmp.*qtmp)./(ntmp.*r.^2)).*wtmp.*zptmp);
        Ih2 = sum((conj(qtmp)./(ntmp.*r.^2)).*wtmp.*zptmp);
    
        sigmac(j) = sigmac(j) - (2*btar(j)*imag(Ic1) + 1i*conj(btar(j)*(Ic2+Ih1-conj(zsrc(j))*Ih2)))/(4*pi);

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

                Ic1 = -sum((qtmp./(ntmp.*r)).*wtmp.*zptmp);
                Ic2 = -sum((conj(qtmp)./(ntmp.*r)).*wtmp.*zptmp);
                Ih1 = sum((conj(ztmp.*qtmp)./(ntmp.*r.^2)).*wtmp.*zptmp);
                Ih2 = sum((conj(qtmp)./(ntmp.*r.^2)).*wtmp.*zptmp);

                sigmac(j) = sigmac(j) - (2*btar(j)*imag(Ic1) + 1i*conj(btar(j)*(Ic2+Ih1-conj(zsrc(j))*Ih2)))/(4*pi);
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

    Ic1 = -cauchy_on_surface_evaluation(qtmp./ntmp,ztmp,zptmp,wtmp,panel_breaks_z,type,periodic_rep);
    Ic2 = -cauchy_on_surface_evaluation(conj(qtmp)./ntmp,ztmp,zptmp,wtmp,panel_breaks_z,type,periodic_rep);
    Ih1 = hypersingular_on_surface_evaluation_zsrc(qtmp./ntmp,ztmp,zptmp,wtmp,panel_breaks_z,type,periodic_rep);
    Ih2 = hypersingular_on_surface_evaluation(conj(qtmp)./ntmp,ztmp,zptmp,wtmp,panel_breaks_z,type,periodic_rep);
    
    sigmac(indices) = sigmac(indices) + (2*btmp.*imag(Ic1) + 1i*conj(btmp.*(Ic2+Ih1-conj(ztmp).*Ih2)))/(4*pi);
    
    wall_start = wall_start + npan + 1;
end