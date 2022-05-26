function Ic = cauchy_on_surface_evaluation(qsrc, zsrc, zpsrc, wsrc, ...
            panel_breaks_z, type)
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
%
%outputs:
% -Ih: the value of the principle value of the integral at the quadrature
% points

sq = special_quad(32);
Nsrc = length(qsrc);
npan = Nsrc/16;

% check if boundary is a "closed curve". Used to handle periodicities.
% Assumes that the panels are of equal length.
panel_length = abs(panel_breaks_z(2)-panel_breaks_z(1));
endpanel_length = abs(panel_breaks_z(end)-panel_breaks_z(1));
if abs(panel_length-endpanel_length) > 1e-12
    closed_curve = 0;
else
    closed_curve = 1;
end

Ic = zeros(Nsrc,1);

for i = 1:Nsrc
    panel_k = ceil(i/16);
    
    switch panel_k
        case 1
            local_panels = [npan, 1, 2];

        case npan - 1
            local_panels = [panel_k - 1, panel_k, panel_k + 1];

        case npan
            local_panels = [npan-1, npan, 1];

        otherwise
            local_panels = [panel_k - 1, panel_k, panel_k + 1];
    end
    
    
    for j = 1:3
        
%         if (j == 1 && local_panels(j) == npan)
%             continue;
%         end
        
        local_indices = 16*(local_panels(j)-1)+1: 16*local_panels(j);
        
        za = panel_breaks_z(local_panels(j));
        if local_panels(j) ~= npan
            zb = panel_breaks_z(local_panels(j)+1);
        else
            if closed_curve
                zb = panel_breaks_z(1);
            else
                % add the difference between the current and previous panel
                zb = panel_breaks_z(local_panels(j)) + panel_breaks_z(local_panels(j))-panel_breaks_z(local_panels(j)-1);
            end
        end
        
        % scale points
        mid = (za + zb)/2;
        len = zb - za;
        
        nzsrc = 2*(zsrc(local_indices) - mid)/len;
        
        if (j == 1 && local_panels(j) == npan)
            mid = panel_breaks_z(1) - (zb - mid);
        end
        
        if (j == 3 && local_panels(j) == 1)
            mid = panel_breaks_z(end) + len + mid - za; %% NOTE: this assumes equal panel length, which might be appropriate
        end
       
        nz = 2*(zsrc(i)-mid)/len;
        p0 = sq.compute_exact_log(nz, nzsrc);

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
    end
    
    % add on regular contribution from points not on adjacent panels
    indices = 1:Nsrc;
    local_indices = [];
    for j = 1:3
        if (j == 1 && local_panels(j) == npan)
            continue;
        end
        
        if (j == 3 && local_panels(j) == 1)
            continue;
        end
        
        local_indices = [local_indices, 16*(local_panels(j)-1)+1: 16*local_panels(j)];
    end
    
    for j = local_indices
        indices(indices == j) = [];
    end
    
    Ic(i) = Ic(i) ...
        + sum(qsrc(indices)./(zsrc(i) - zsrc(indices)).*zpsrc(indices).*wsrc(indices));
end
