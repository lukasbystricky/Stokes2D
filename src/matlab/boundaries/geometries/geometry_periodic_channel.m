function [z,zp,zpp] = geometry_periodic_channel(y,yp,ypp,T,dir,Lx)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Function to compute a periodic channel wall. The normals of the wall 
% must be pointing into the channel, so the direction is important.
%
% Input:
%       y, function for the position 
%       yp, function for the first derivative
%       ypp, function for the second derivative
%       T, parametrisation [0,2*pi]
%       dir, which direction the parametrisation goes (1 for top wall, -1
%               for bottom wall)
%       Lx, periodic length scale in the x direction
%
% Output:
%       z, function returning position at t
%       zp, function returning first derivative at t
%       zpp, function returning second derivative at t
%       
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% map T to interval [0, L]
s = T * Lx /(2*pi);

if dir > 0
    z = s-(Lx/2)+1i*y(s);
    zp = Lx/(2*pi)*(1+1i*yp(s));
    zpp = (Lx/(2*pi))^2*1i*ypp(s);
else
    z = (Lx-s)-Lx/2+1i*y(Lx-s);
    zp = -(Lx/(2*pi))*(1+1i*yp(Lx-s));
    zpp = (Lx/(2*pi))^2*1i*ypp(Lx-s);
end