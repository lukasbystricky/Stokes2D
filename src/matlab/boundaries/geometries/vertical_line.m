function [z,zp,zpp] = vertical_line(T,x,ymin,ymax,dir)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Function to compute a simple vertical line. Used in the pouseille flow
% averages computations. 
%
% Input:
%       T, parametrisation [0,2*pi]
%       x, position along the x axis
%       ymin, smallest value along the y axis
%       ymax, largest value along the y axis
%       dir, which direction the parametrisation goes (1 from ymin to ymax
%           and -1 from ymax to ymin)
%
% Output:
%       z, function returning position at t
%       zp, function returning first derivative at t
%       zpp, function returning second derivative at t
%       
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
L = ymax-ymin;
if L < 0
    error('ymax must not be less than ymin');
end
if min(T) < 0 || max(T) > 2*pi
    error('T must be in [0,2*pi]');
end
s = T * L /(2*pi);

if dir > 0
    z = x + 1i*(s-L/2);
    zp = L/(2*pi)*1i;
else
    z = x - 1i*(s-L/2);
    zp = -L/(2*pi)*1i;
end
zpp = 0;