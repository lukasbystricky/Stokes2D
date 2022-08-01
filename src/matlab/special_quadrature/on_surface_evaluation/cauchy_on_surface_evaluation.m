function Ic = cauchy_on_surface_evaluation(qsrc, zsrc, zpsrc, wsrc, ...
            panel_breaks_z, type, periodic_rep)
%CAUCHY_ON_SURFACE_EVALUATION evaluates the principal value of the
%integral q(tau)/(z_i - tau), where z_i coincides with the quadrature 
%points on the boundary. Corrects the value using special quadrature for
%points on the same panel and adjacent panels as the target point. Includes
%the limiting value as the target point approaches the boundary from the
%outside. 
%
%inputs:
% -qsrc: density at quadrature points
% -zsrc: quadrature points
% -zpsrc: dz/dt at quadrature points
% -wsrc: quadrature weights
% -panel_breaks_z: panel endpoints in physical space
% -periodic_rep: boolean indicating if periodic replicates should be 
% compensated for
%
%outputs:
% -Ic: the value of the principle value of the integral at the quadrature
% points

sq = special_quad(32);
Nsrc = length(qsrc);
npan = Nsrc/16;

Ic = zeros(Nsrc,1);

for i = 1:Nsrc
    panel_k = ceil(i/16);
    
    switch panel_k
        case 1
            local_panels = [npan, 1, 2];
            next_adjacent_panel = 3;

        case npan - 1
            local_panels = [panel_k - 1, panel_k, panel_k + 1];
            next_adjacent_panel = 1;

        case npan
            local_panels = [npan-1, npan, 1];
            next_adjacent_panel = 2;

        otherwise
            local_panels = [panel_k - 1, panel_k, panel_k + 1];
            next_adjacent_panel = panel_k + 2;
    end
    
    
    for j = 1:3
        
        local_indices = 16*(local_panels(j)-1)+1: 16*local_panels(j);
        mask = 1:16;
        
        % endpoints of panel
        za = panel_breaks_z(local_panels(j));
        zb = panel_breaks_z(local_panels(j)+1);

        % scale points
        mid = (za + zb)/2;
        len = zb - za;
        
        nzsrc = 2*(zsrc(local_indices) - mid)/len;
        
        if periodic_rep
            if (j == 1 && local_panels(j) == npan)
                mid = panel_breaks_z(1) - (zb - mid);
            end
            if (j == 3 && local_panels(j) == 1)
                mid = panel_breaks_z(end) + mid - za;
            end
            
            % correct only the closest four points on each adjacent panel
%             if j == 1
%                 mask = 13:16;
%             end
%             if j == 3
%                 mask = 1:4;
%             end
        end
        
        nz = 2*(zsrc(i)-mid)/len;
        p0 = sq.compute_exact_log(nz, nzsrc);

%         test_ind = 16*(next_adjacent_panel-1)+1: 16*next_adjacent_panel;
%         test_za = panel_breaks_z(next_adjacent_panel);
%         test_zb = panel_breaks_z(next_adjacent_panel+1);
%         exact_int = log(test_zb-zsrc(i)) - log(test_za-zsrc(i));
%         [accurate, testsum, err] = sq.sq_necessary(exact_int, zsrc(i), zsrc(test_ind), ...
%             zpsrc(test_ind), wsrc(test_ind));
        %assert(accurate,'%e',err);

        if j == 2 
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
        Ic(i) = Ic(i) ...
            +(sq.cauchy_integral(qsrc(local_indices), nz, nzsrc, p0));
%         Ic(i) = Ic(i) ...
%             +(sq.cauchy_integral_16(qsrc(local_indices), nz, nzsrc, p0, mask));
    end
    
    % add on regular contribution from points not on adjacent panels
    indices = 1:Nsrc;
    local_indices = [];
    for j = 1:3
        if periodic_rep
            if (j == 1 && local_panels(j) == npan)
                continue;
            end

            if (j == 3 && local_panels(j) == 1)
                continue;
            end
        end
        
        tmp = 16*(local_panels(j)-1)+1: 16*local_panels(j);
        tmp_ind = tmp;
        
        % correct only the closest four points on each adjacent panel
%         if j == 1
%             tmp_ind = tmp(end-3:end);
%         end
%         if j == 3
%             tmp_ind = tmp(1:4);
%         end

        local_indices = [local_indices, tmp_ind];
    end
    
    for j = local_indices
        indices(indices == j) = [];
    end

    Ic(i) = Ic(i) ...
        + sum(qsrc(indices)./(zsrc(i) - zsrc(indices)).*zpsrc(indices).*wsrc(indices));
end
