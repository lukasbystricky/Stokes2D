function ztar = find_target_pt(z,mid,bnds)
%FIND_TARGET_PT finds the periodic replicate of the target point closest to
%the panel.
%
%inputs:
% -z: target point in complex variables
% -mid: middle point of source panel
% -bnds: bounds of periodic box (reference cell)
%
%outputs:
% -ztar: the periodic replicate of z that is closest to the panel

ztar = z;

d0 = abs(ztar-mid);
d1 = abs(mid - (ztar - (bnds(2)-bnds(1))));
d2 = abs(mid - (ztar + (bnds(2)-bnds(1))));
d3 = abs(mid - (ztar - 1i*(bnds(4)-bnds(3))));
d4 = abs(mid - (ztar + 1i*(bnds(4)-bnds(3))));
d5 = abs(mid - (ztar + ((bnds(2)-bnds(1))+1i*(bnds(4)-bnds(3)))));
d6 = abs(mid - (ztar + (-(bnds(2)-bnds(1))+1i*(bnds(4)-bnds(3)))));
d7 = abs(mid - (ztar + ((bnds(2)-bnds(1))-1i*(bnds(4)-bnds(3)))));
d8 = abs(mid - (ztar - ((bnds(2)-bnds(1))+1i*(bnds(4)-bnds(3)))));

[~,idx] = min([d0 d1 d2 d3 d4 d5 d6 d7 d8]);
periodic_case = idx - 1;

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
