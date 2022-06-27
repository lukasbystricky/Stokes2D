function I = on_surface_evaluation(p,qsrc,zsrc,zpsrc,wsrc,panel_breaks_z,...
    lim_dir,nbr_neighbor_pts,periodic_rep,varargin)
%ON_SURFACE_EVALUATION evaluates the principal value of the
%integral q(tau)/(z_i - tau)^p, where z_i coincides with the quadrature 
%points on the boundary. Corrects the value using special quadrature for
%points on the same panel and adjacent panels as the target point. Includes
%the limiting value as the target point approaches the boundary from the
%outside. 
%
%inputs:
% -p: power of denominator (1, 2, or 3)
% -qsrc: density at quadrature points
% -zsrc: quadrature points
% -zpsrc: dz/dt at quadrature points
% -wsrc: quadrature weights
% -panel_breaks_z: panel endpoints in physical space
% -lim_dir: direction of the limiting value
% -nbr_neighbor_pts: number of neighboring points to correct
% -periodic_rep: boolean indicating if periodic replicates should be 
% compensated for
%
%outputs:
% -I: the principle value of the integral at the quadrature points

if nargin > 9
    reference_cell = varargin{1};
end

% check validity of input
assert(length(qsrc) == length(zsrc),'dimension of arrays not consistent');
assert(p == 1 || p == 2 || p == 3,'incorrect p: must be 1, 2, or 3');
assert(any(strcmp(lim_dir,{'fluid','solid','surface'})),...
    'incorrect lim_dir: must be "fluid", "solid", or "surface"');
assert(periodic_rep && (nargin > 9),'reference cell is missing');

% setup
sq = special_quad(32);
Nsrc = length(qsrc);
panels_per_wall = Nsrc/16;
I = zeros(Nsrc,1);

% loop over each panel on current wall
for j = 1:panels_per_wall
    
    % source and target indices
    source_ind = (j-1)*16+(1:16);
    target_ind = mod((j-1)*16+(-(nbr_neighbor_pts-1):nbr_neighbor_pts+16)+panels_per_wall*16-1, panels_per_wall*16) + 1;
    
    % endpoints of panel
    za = panel_breaks_z(j);
    zb = panel_breaks_z(j+1);

    % middle and length of panel
    mid = (za + zb)/2;
    len = zb - za;

    % scaled source panel
    nzsrc = 2*(zsrc(source_ind)-mid)/len;

    % quantities on source panel
    zsrc_tmp = zsrc(source_ind);
    zpsrc_tmp = zpsrc(source_ind);
    wsrc_tmp = wsrc(source_ind);
    qsrc_tmp = qsrc(source_ind);

    % loop over each target point
    for k = 1:length(target_ind)
        
        % local target index and what panel it lies on
        tar_ind_tmp = target_ind(k);
        tar_pan = ceil(tar_ind_tmp/16);
        
        % find correct target point
        if periodic_rep
            [ztar_tmp,~] = find_target_pt(zsrc(tar_ind_tmp),mid,len,reference_cell);
        else
            ztar_tmp = zsrc(tar_ind_tmp);
        end
        
        % scaled target point
        nz = 2*(ztar_tmp-mid)/len;
        
        % initiate recursion
        p0 = sq.compute_exact_log(nz, nzsrc);

        % add limiting value if needed
        if j == tar_pan
            switch lim_dir
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
        
        % apply product integration correction (inconsistent signs due to
        % the implementation in special_quad.m)
        if p == 1
            I(tar_ind_tmp) = I(tar_ind_tmp) - sq.cauchy_integral(qsrc_tmp,nz,nzsrc,p0);
        elseif p == 2
            I(tar_ind_tmp) = I(tar_ind_tmp) + sq.hypersingular_integral(qsrc_tmp,nz,nzsrc,p0,zb,za);
        elseif p == 3
            I(tar_ind_tmp) = I(tar_ind_tmp) - sq.supersingular_integral(qsrc_tmp,nz,nzsrc,p0zb,za);
        end
        
    end
    
    % add on regular contribution from non-neighboring points
    target_far_ind = setdiff(1:Nsrc,target_ind);
    
    if j == 1
        target_far_ind = [target_far_ind (Nsrc-(nbr_neighbor_pts-1):Nsrc)];
    end
    if j == panels_per_wall
        target_far_ind = [1:nbr_neighbor_pts target_far_ind];
    end
    
    for k = target_far_ind
        I(k) = I(k) - sum(qsrc_tmp./(zsrc(k) - zsrc_tmp).^p.*zpsrc_tmp.*wsrc_tmp);
    end
end
