function [z,zp,zpp] = geometry_star(A,k,r,center,angle,T,varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% A boundary function for star-shaped particle.
%
% Input:
%       A, magnitude of star legs
%       k, number of star legs
%       r, radius of a leg from the particle center
%       center, center of star
%       angle, rotation of the star in positive direction 
%       T, parametrisation [0,2*pi]
%
% Output:
%       z, function returning position at t
%       zp, function returning first derivative at t
%       zpp, function returning second derivative at t
%       
% Updated: - Anna Jan 10, 2023
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if nargin > 6
    if varargin{1} == 1
        dir = 1;
    else 
        dir = -1;
    end
else
    dir = 1;
end

s = r/(1+A);    % set scaling of the particle to match certain radius

fz = @(t) (1+A*cos(k*t)).*exp(1i*t);
fzz = @(t) (-A*k*sin(k*t)).*exp(1i*t)+1i*fz(t);
fzzz = @(t) exp(1i*(t)).*(-1+(-A*k^2-A)*cos(k*(t))-2*A*k*1i*sin(k*(t)));

% NB: think about the direction of the rotation angle here
if dir == 1
    z = s*fz(T+angle)+center;   % scale star
    zp = s*fzz(T+angle);
    zpp = s*fzzz(T+angle);
else
    z = s*fz(-T+angle)+center;
    zp = -s*fzz(-T+angle);
    zpp = s*fzzz(-T+angle);
end
