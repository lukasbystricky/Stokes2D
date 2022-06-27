function Pc = pressure_slp_on_surface_correction_new(P, solution_local, lim_dir)
%PRESSURE_SLP_ON_SURFACE_CORRECTION corrects the single-layer 
%pressure for on-surface evaluation. Adds on the jump as the target
%approaches the boundary from different directions.
%
%inputs:
% -P: uncorrected single-layer pressure evaluated at quadrature points on 
% surface
% -solution:local: solution structure containing geometry information, as 
% well as the density 
% -lim_dir: string describing from what direction the target is approaching
% the boundary
%
%output:
% -Pc: corrected pressure

% geometry information of local solution
domain_local = solution_local.problem.domain;
periodic = solution_local.problem.periodic;
reference_cell = domain_local.reference_cell;
nbr_neighbor_pts = domain_local.nbr_neighboring_points;

% quantities from local solution
zsrc = domain_local.z;
zpsrc = domain_local.zp;
nsrc = -1i*zpsrc./abs(zpsrc);
wsrc = domain_local.quad_weights;
qsrc = solution_local.q(:,1) + 1i*solution_local.q(:,2);

% initialize indices
wall_indices = domain_local.wall_indices;
wall_start = 1;

% initialize corrected quantity
Pc = P;

% loop over every wall
for i = 1:size(domain_local.wall_indices,1)
    
    % information about current wall
    thiswall = wall_indices(i,1):wall_indices(i,2);
    panels_per_wall = (wall_indices(i,2)-wall_indices(i,1)+1)/16;
    panel_breaks_z = domain_local.panel_breaks(wall_start:wall_start + panels_per_wall);

    % check if boundary is a closed curve
    if abs(panel_breaks_z(1)-panel_breaks_z(end)) < 1e-12
        closed_curve = 1;
    else
        closed_curve = 0;
    end
    
    % if problem is periodic and the domain is not a closed curve, then we
    % need to handle periodic replicates in a specific way
    periodic_rep = periodic && ~closed_curve;

    % subract off contribution from the current wall
    for j = thiswall

        indices_tmp = thiswall;
        indices_tmp(thiswall==j) = [];
        
        ztmp = zsrc(indices_tmp);
        qtmp = qsrc(indices_tmp);
        wtmp = wsrc(indices_tmp);
        ntmp = nsrc(indices_tmp);
        zptmp = zpsrc(indices_tmp);
        r = zsrc(j) - ztmp;
        
        Ic = -sum((qtmp./(ntmp.*r)).*wtmp.*zptmp);
        
        Pc(j) = Pc(j) - -imag(Ic)/(2*pi);

    end
    
    % subtract of contribution from periodic replicates due to them being
    % included in Ewald, will be computed later accurately using a special
    % quadrature rule
    if periodic_rep
        for j = [1 panels_per_wall]
                        
            % source and target indices
            source_ind = wall_indices(i,1)-1 + ((j-1)*16+(1:16));
            target_ind = wall_indices(i,1) + ...
                mod((j-1)*16+(-(nbr_neighbor_pts-1):nbr_neighbor_pts+16)+panels_per_wall*16-1, panels_per_wall*16);
            
            % consider only the closest nbr_neighbor_pts target replica pts
            if j == 1
                target_ind = target_ind(1:nbr_neighbor_pts);
            else
                target_ind = target_ind(end-(nbr_neighbor_pts-1):end);
            end
            
            % quantities on source panel
            ztmp = zsrc(source_ind);
            qtmp = qsrc(source_ind);
            wtmp = wsrc(source_ind);
            ntmp = nsrc(source_ind);
            zptmp = zpsrc(source_ind);
            
            % endpoints of panel
            za = panel_breaks_z(j);
            zb = panel_breaks_z(j+1);

            % middle and length of panel
            mid = (za + zb)/2;
            len = zb - za;
            
            for k = 1:length(target_ind)
                
                % find periodic target point
                tar_ind_tmp = target_ind(k);
                [ztar_tmp,~] = find_target_pt(zsrc(tar_ind_tmp),mid,len,reference_cell);
                r = ztar_tmp - ztmp;

                Ic = -sum((qtmp./(ntmp.*r)).*wtmp.*zptmp);
        
                Pc(tar_ind_tmp) = Pc(tar_ind_tmp) - -imag(Ic)/(2*pi);
                
            end
        end
    end

    % evaluate integral using special quadrature
    Ic = on_surface_evaluation(1,qsrc(thiswall)./nsrc(thiswall),...
        zsrc(thiswall),zpsrc(thiswall),wsrc(thiswall),panel_breaks_z,...
        lim_dir,nbr_neighbor_pts,periodic_rep,reference_cell);

    % add on correction
    Pc(thiswall) = Pc(thiswall) + -imag(Ic)/(2*pi); 

    % update panel count
    wall_start = wall_start + panels_per_wall + 1;
    
end
