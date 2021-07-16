function [z,zp,zpp,zppp] = flower_func(t)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Evaluates a Narms-petaled flower with petal length 0.5*a.
%
% Input:
%       a, petal length*2 
%       scale, scaling term for size of flower
%       theta, tilt of flower
%       Narms, number of petals of flowers
%       center, center for flower
%       phi, parametrisation [0,2*pi]
%
% Output:
%       z, function returning position at t
%       zp, function returning first derivative at t
%       zpp, function returning second derivative at t
%       
% Updated: -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

s = 0;
a = 0.3;
b = 5;
z = (1+a*cos(b*(t+s))).*exp(1i*(t+s)); %starfish parametrization
zp = (-a*b*sin(b*(t+s))+1i*(1+a*cos(b*(t+s)))).*exp(1i*(t+s));
zpp = exp(1i*(t+s)).*(-1+(-a*b^2-a)*cos(b*(t+s))-2*a*b*1i*sin(b*(t+s)));
zppp = exp(1i*(t+s)).*(1i*(-1+(-a*b^2-a)*cos(b*(t+s))-2*a*b*1i*sin(b*(t+s))) + ...
    -b*(-a*b^2-a)*sin(b*(t+s)) - 2*a*b^2*1i*cos(b*(t+s)));
        
