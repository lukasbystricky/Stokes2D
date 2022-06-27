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
periodic = solution_local.problem.periodic;

zsrc = domain.z;
zpsrc = domain.zp;
nsrc = -1i*zpsrc./abs(zpsrc);
wsrc = domain.quad_weights;

qsrc = (solution_local.q(:,1) + 1i*solution_local.q(:,2))./nsrc;

Pc = P;
Ic_old = [];
Pc_old = zeros(length(zsrc),1);

wall_start = 1;
for i = 1:size(domain.wall_indices,1)

    indices = domain.wall_indices(i,1):domain.wall_indices(i,2);
    npan = length(indices)/16;

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
        
        Pc(j) = Pc(j) - imag(sum(qsrc(indices_tmp).*wsrc(indices_tmp).*...
                zpsrc(indices_tmp)./(zsrc(j) - zsrc(indices_tmp))))/(2*pi);

        if periodic_rep
            % we also correct for the adjacent panels that have been
            % periodically replicated later, so these need to be removed too
            indices_tmp = [];

            if (j - (i-1)*16*npan) <= 16
                indices_tmp = indices(end-15:end);
                %indices_tmp = indices(end-3:end);
                
                if real(zsrc(indices(2))) > real(zsrc(indices(1)))
                    ztmp = zsrc(indices_tmp) - solution_local.problem.Lx;
                else
                    ztmp = zsrc(indices_tmp) + solution_local.problem.Lx;
                end

            elseif (j - (i-1)*16*npan) >= (length(indices) - 15)
                indices_tmp = indices(1:16);
                %indices_tmp = indices(1:4);

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
                zptmp = zpsrc(indices_tmp);

                Ic = -sum((qtmp./r).*wtmp.*zptmp);
                
                Pc(j) = Pc(j) - -imag(Ic)/(2*pi);

            end
        end
    end    
    Pc_old(indices) = Pc(indices);
    %add on special quadrature
    Ic = -cauchy_on_surface_evaluation(qsrc(indices), zsrc(indices), zpsrc(indices), wsrc(indices), panel_breaks_z, type, periodic_rep);
    Pc(indices) = Pc(indices) + -imag(Ic)/(2*pi); 
    Ic_old = [Ic_old; Ic];
    wall_start = wall_start + npan + 1;       
end

%%%%
domain_local = solution_local.problem.domain;
Pc = P;
sq = special_quad(32);
Nsrc = length(zsrc);
wall_indices = domain_local.wall_indices;
wall_start = 1;
Ic = zeros(Nsrc,1);

for i = 1:size(domain_local.wall_indices,1)

    thiswall = wall_indices(i,1):wall_indices(i,2);
    panels_per_wall = (wall_indices(i,2)-wall_indices(i,1)+1)/16;
    panel_breaks_z = domain.panel_breaks(wall_start:wall_start + panels_per_wall);
    
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
    for j = thiswall
       
        indices_tmp = thiswall;
        indices_tmp(thiswall==j) = [];
        
        Pc(j) = Pc(j) - imag(sum(qsrc(indices_tmp).*wsrc(indices_tmp).*...
                zpsrc(indices_tmp)./(zsrc(j) - zsrc(indices_tmp))))/(2*pi);

    end
    
    if periodic_rep
        for j = [1 panels_per_wall]

            % source and target indices
            source_ind = wall_indices(i,1)-1+((j-1)*16+(1:16));
            target_ind = wall_indices(i,1)+...
                mod((j-1)*16+(-3:20)+panels_per_wall*16-1, panels_per_wall*16);
            
            if j == 1
                target_ind = target_ind(1:4);
            else
                target_ind = target_ind(end-3:end);
            end
        
            zsrctmp = zsrc(source_ind);
            qtmp = qsrc(source_ind);
            wtmp = wsrc(source_ind);
            zptmp = zpsrc(source_ind);
            
            for k = 1:length(target_ind)

                tar_ind_tmp = target_ind(k);
                tar_pan = ceil(tar_ind_tmp/16-(i-1)*panels_per_wall);
                ztar_tmp = zsrc(tar_ind_tmp);

                if j == 1 && tar_pan == panels_per_wall
                    if real(zsrc(thiswall(2))) > real(zsrc(thiswall(1)))
                        ztar_tmp = zsrc(tar_ind_tmp) - solution_local.problem.Lx;
                    else
                        ztar_tmp = zsrc(tar_ind_tmp) + solution_local.problem.Lx;
                    end
                end

                if j == panels_per_wall && tar_pan == 1
                    if real(zsrc(thiswall(2))) > real(zsrc(thiswall(1)))
                        ztar_tmp = zsrc(tar_ind_tmp) + solution_local.problem.Lx;
                    else
                        ztar_tmp = zsrc(tar_ind_tmp) - solution_local.problem.Lx;
                    end
                end
                
                Pc(tar_ind_tmp) = Pc(tar_ind_tmp) - imag(sum(qtmp.*wtmp.*...
                    zptmp./(ztar_tmp - zsrctmp)))/(2*pi);
            
            end
        end
    end
        
    for j = 1:panels_per_wall
        
        % source and target indices
        source_ind = wall_indices(i,1)-1+((j-1)*16+(1:16));
        target_ind = wall_indices(i,1)+...
            mod((j-1)*16+(-3:20)+panels_per_wall*16-1, panels_per_wall*16);
        
        % endpoints of panel
        za = panel_breaks_z(j);
        zb = panel_breaks_z(j+1);

        % scale points
        mid = (za + zb)/2;
        len = zb - za;

        % scaled source panel
        nzsrc = 2*(zsrc(source_ind)-mid)/len;

        % q source panel
        zsrc_tmp = zsrc(source_ind);
        zpsrc_tmp = zpsrc(source_ind);
        wsrc_tmp = wsrc(source_ind);
        qsrc_tmp = qsrc(source_ind);

        % loop over each target point
        for k = 1:length(target_ind)
            
            tar_ind_tmp = target_ind(k);
            tar_pan = ceil(tar_ind_tmp/16-(i-1)*panels_per_wall);
            ztar_tmp = zsrc(tar_ind_tmp);
            
            if j == 1 && tar_pan == panels_per_wall
                if real(zsrc(thiswall(2))) > real(zsrc(thiswall(1)))
                    ztar_tmp = zsrc(tar_ind_tmp) - solution_local.problem.Lx;
                else
                    ztar_tmp = zsrc(tar_ind_tmp) + solution_local.problem.Lx;
                end
            end
            
            if j == panels_per_wall && tar_pan == 1
                if real(zsrc(thiswall(2))) > real(zsrc(thiswall(1)))
                    ztar_tmp = zsrc(tar_ind_tmp) + solution_local.problem.Lx;
                else
                    ztar_tmp = zsrc(tar_ind_tmp) - solution_local.problem.Lx;
                end
            end

            nz = 2*(ztar_tmp-mid)/len;
            p0 = sq.compute_exact_log(nz, nzsrc);
            
            %if k > 4 && k < 17
            if j == tar_pan
                switch type % add limiting value if needed
                    case 'fluid'
                        if imag(nz) > 0
                            p0 = p0 - 2*1i*pi;
                        end

                    case 'solid'
                        if imag(nz) < 0
                            p0 = p0 + 2*1i*pi;
                        end

                    case 'surface'
                        if imag(nz) > 0
                            p0 = p0 - 1i*pi;
                        else
                            p0 = p0 + 1i*pi;
                        end
                end
            end
            
            % apply product integration correction
            Ic(tar_ind_tmp) = Ic(tar_ind_tmp) - sq.cauchy_integral(qsrc_tmp, nz, nzsrc, p0);

        end
        
        % add on regular contribution from all other target points
        target_far_ind = setdiff(domain.wall_indices(i,1):domain.wall_indices(i,2),target_ind);
        
        if j == 1
            target_far_ind = [target_far_ind (domain.wall_indices(i,2)-3:domain.wall_indices(i,2))];
        end
        if j == panels_per_wall
            target_far_ind = [domain.wall_indices(i,1):domain.wall_indices(i,1)+3 target_far_ind];
        end
        
        for k = target_far_ind
            Ic(k) = Ic(k) ...
                - sum(qsrc_tmp./(zsrc(k) - zsrc_tmp).*zpsrc_tmp.*wsrc_tmp);
        end
    end

    wall_start = wall_start + panels_per_wall + 1;
    
end

correction = -imag(Ic)/(2*pi);
Pc = Pc + correction;
