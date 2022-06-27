function [ztar,check_sq] = find_target_pt(z,mid,len,bnds)
%FIND_TARGET_PT finds the periodic replicate of the target point closest to
%the panel. Periodic case -1 means no need for special quadrature, 0 means
%z is close to a panel in the reference cell.
%
%inputs:
% -z: target point in complex variables
% -mid: middle point of source panel
% -len: length of source panel
% -bnds: bounds of periodic box (reference cell)
%
%outputs:
% -ztar: the periodic replicate of z that is closest to the panel
% check_sq: boolean indicating if special quadrature is needed

periodic_case = -1;
ztar = z;

if abs(ztar-mid) < abs(len)
    periodic_case = 0;
elseif abs(mid - (ztar - (bnds(2)-bnds(1)))) < abs(len)
    periodic_case = 1;
elseif abs(mid - (ztar + (bnds(2)-bnds(1)))) < abs(len)
    periodic_case = 2;
elseif abs(mid - (ztar - 1i*(bnds(4)-bnds(3)))) < abs(len)
    periodic_case = 3;
elseif abs(mid - (ztar + 1i*(bnds(4)-bnds(3)))) < abs(len)
    periodic_case = 4;
elseif abs(mid - (ztar + ((bnds(2)-bnds(1))+1i*(bnds(4)-bnds(3))))) < abs(len)
    periodic_case = 5;
elseif abs(mid - (ztar + (-(bnds(2)-bnds(1))+1i*(bnds(4)-bnds(3))))) < abs(len)
    periodic_case = 6;
elseif abs(mid - (ztar + ((bnds(2)-bnds(1))-1i*(bnds(4)-bnds(3))))) < abs(len)
    periodic_case = 7;
elseif abs(mid - (ztar - ((bnds(2)-bnds(1))+1i*(bnds(4)-bnds(3))))) < abs(len)
    periodic_case = 8;
end

switch periodic_case
    case 1
        ztar = ztar - (bnds(2)-bnds(1));
    case 2
        ztar = ztar + (bnds(2)-bnds(1));
    case 3
        ztar = ztar - 1i*(bnds(4)-bnds(3));
    case 4
        ztar = ztar + 1i*(bnds(4)-bnds(3));
    case 5
        ztar = ztar + ((bnds(2)-bnds(1))+1i*(bnds(4)-bnds(3)));
    case 6
        ztar = ztar + (-(bnds(2)-bnds(1))+1i*(bnds(4)-bnds(3)));
    case 7
        ztar = ztar + ((bnds(2)-bnds(1))-1i*(bnds(4)-bnds(3)));
    case 8
        ztar = ztar - ((bnds(2)-bnds(1))+1i*(bnds(4)-bnds(3)));
end

check_sq = periodic_case > -1;
