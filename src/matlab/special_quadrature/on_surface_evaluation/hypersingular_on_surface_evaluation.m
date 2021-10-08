function Ih = hypersingular_on_surface_evaluation(qsrc, zsrc, zpsrc, wsrc,...
                                panel_breaks_z, type)
%HYPERSINGULAR_ON_SURFACE_EVALUATION evaluates the finite part of the
%integral q(tau)/(z_i - tau)^2, where z_i coincides with the quadrature 
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
% -Ih: the value of the finite part integral at the quadrature points

sq = special_quad(32);
Nsrc = length(qsrc);
npan = Nsrc/16;

Ih = zeros(Nsrc,1);

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
        local_indices = 16*(local_panels(j)-1)+1: 16*local_panels(j);
        
        za = panel_breaks_z(local_panels(j));
        if local_panels(j) ~= npan
            zb = panel_breaks_z(local_panels(j)+1);
        else
            zb = panel_breaks_z(local_panels(j)) + (panel_breaks_z(local_panels(j))-panel_breaks_z(local_panels(j)-1));
        end
        
        % scale points
        mid = (za + zb)/2;
        len = zb - za;
        nzsrc = 2*(zsrc(local_indices) - mid)/len;
        
        nz = 2*(zsrc(i)-mid)/len;
        
        p0 = sq.compute_exact_log(nz, nzsrc);
        
        % apply product integration correction
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
        
        Ih(i) = Ih(i) ...
            + sq.hypersingular_integral(qsrc(local_indices), nz, nzsrc, ...
                p0, zb, za);        
    end
    
    indices = 1:Nsrc;
    local_indices = [];
    for j = 1:3
        local_indices = [local_indices, 16*(local_panels(j)-1)+1: 16*local_panels(j)];
    end
    
    for j = local_indices
        indices(indices == j) = [];
    end
    
    Ih(i) = Ih(i) ...
        + sum(qsrc(indices)./(zsrc(indices) - zsrc(i)).^2.*zpsrc(indices).*wsrc(indices));
end
