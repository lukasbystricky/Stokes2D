function dPc = pressure_gradient_dlp_on_surface_correction_new(dP,solution_local,lim_dir)
%PRESSURE_GRADIENT_DLP_ON_SURFACE_CORRECTION corrects the double-layer 
%pressure gradient for on-surface evaluation. Adds on the jump as the target
%approaches the boundary from different directions.
%
%inputs:
% -dP: uncorrected double-layer pressure gradient evaluated at quadrature points on 
% surface
% -solution:local: solution structure containing geometry information, as 
% well as the density 
% -lim_dir: string describing from what direction the target is approaching
% the boundary
%
%output:
% -dPc: corrected pressure gradient

% geometry information of local solution
domain_local = solution_local.problem.domain;
periodic = solution_local.problem.periodic;
reference_cell = domain_local.reference_cell;
nbr_neighbor_pts = domain_local.nbr_neighbor_pts;

% quantities from local solution
zsrc = domain_local.z;
zpsrc = domain_local.zp;
wsrc = domain_local.quad_weights;
qsrc = solution_local.q(:,1) + 1i*solution_local.q(:,2);

% initialize indices
wall_indices = domain_local.wall_indices;
wall_start = 1;

% initialize corrected quantity
dPc = dP;

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
        zptmp = zpsrc(indices_tmp);
        r = zsrc(j) - ztmp;
        
        Is = -sum((qtmp./r.^3).*wtmp.*zptmp);
        
        dPc(j) = dPc(j) - 1i*conj(Is)/pi;

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
            zptmp = zpsrc(source_ind);
            
            % endpoints of panel
            za = panel_breaks_z(j);
            zb = panel_breaks_z(j+1);

            % middle of panel
            mid = (za + zb)/2;
            
            for k = 1:length(target_ind)
                
                % find periodic target point
                tar_ind_tmp = target_ind(k);
                ztar_tmp = find_target_pt(zsrc(tar_ind_tmp),mid,reference_cell);
                r = ztar_tmp - ztmp;

                Is = -sum((qtmp./r.^3).*wtmp.*zptmp);
        
                dPc(tar_ind_tmp) = dPc(tar_ind_tmp) - 1i*conj(Is)/pi;
                
            end
        end
    end

    % quantities on this wall
    ztmp = zsrc(thiswall);
    qtmp = qsrc(thiswall);
    wtmp = wsrc(thiswall);
    zptmp = zpsrc(thiswall);
    
    % evaluate integral using special quadrature
    Is = on_surface_evaluation(3,qtmp,ztmp,zptmp,wtmp,panel_breaks_z,...
        lim_dir,nbr_neighbor_pts,periodic_rep,reference_cell);

    % add on correction
    dPc(thiswall) = dPc(thiswall) + 1i*conj(Is)/pi;

    % update panel count
    wall_start = wall_start + panels_per_wall + 1;
    
end
