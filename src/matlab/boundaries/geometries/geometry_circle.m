function [z,zp,zpp] = geometry_circle(R,center,T,varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% A boundary function for a circle.
%
% Input:
%       R, radius of circle 
%       center, center of circle
%       T, parametrisation [0,2*pi]
%
% Output:
%       z, function returning position at t
%       zp, function returning first derivative at t
%       zpp, function returning second derivative at t
%       
% Updated: - David Jan 10, 2023
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if nargin > 3
    if varargin{1} == 1
        dir = 1;
    else 
        dir = -1;
    end
else
    dir = 1;
end

if dir == 1
    z = R*exp(1i*T)+center;
    zp = 1i*R*exp(1i*T);
    zpp = -R*exp(1i*T);
else
    z = R*exp(-1i*T)+center;
    zp = -1i*R*exp(-1i*T);
    zpp = -R*exp(-1i*T);
end   