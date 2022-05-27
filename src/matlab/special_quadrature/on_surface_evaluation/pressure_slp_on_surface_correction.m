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
            % We will also correct for the adjacent panels that have been
            % periodically replicated later, so these need to be removed too

            if (j - (i-1)*16*npan <= 16)
               indices_tmp = indices(end-15:end);

               if (real(zsrc(indices(2))) > real(zsrc(indices(1))))
                   Pc(j) = Pc(j) - imag(sum(qsrc(indices_tmp).*wsrc(indices_tmp).*...
                    zpsrc(indices_tmp)./(zsrc(j) - (zsrc(indices_tmp) - solution_local.problem.Lx))))/(2*pi); %% NOTE: this only works for walls that are periodic in the x-direction!
               else
                   Pc(j) = Pc(j) - imag(sum(qsrc(indices_tmp).*wsrc(indices_tmp).*...
                    zpsrc(indices_tmp)./(zsrc(j) - (zsrc(indices_tmp) + solution_local.problem.Lx))))/(2*pi);
               end
            end

            if (j - (i-1)*16*npan >= length(indices) - 15)
               indices_tmp = indices(1:16);

               if (real(zsrc(indices(2))) > real(zsrc(indices(1))))
                  Pc(j) = Pc(j) - imag(sum(qsrc(indices_tmp).*wsrc(indices_tmp).*...
                    zpsrc(indices_tmp)./(zsrc(j) - (zsrc(indices_tmp) + solution_local.problem.Lx))))/(2*pi); %% NOTE: this only works for walls that are periodic in the x-direction!
               else
                  Pc(j) = Pc(j) - imag(sum(qsrc(indices_tmp).*wsrc(indices_tmp).*...
                    zpsrc(indices_tmp)./(zsrc(j) - (zsrc(indices_tmp) - solution_local.problem.Lx))))/(2*pi); 
               end
            end
        end

    end    
    
    %add on special quadrature
    Ic = cauchy_on_surface_evaluation(qsrc(indices), zsrc(indices), zpsrc(indices), wsrc(indices), panel_breaks_z, type, periodic_rep);
    Pc(indices) = Pc(indices) + imag(Ic)/(2*pi); 
            
    wall_start = wall_start + npan + 1;       
end