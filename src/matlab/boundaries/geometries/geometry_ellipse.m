function [z,zp,zpp] = geometry_ellipse(a,b,center,theta,T,varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% A boundary function for an ellipse.
%
% Input:
%       a, major semi-axis of ellipse
%       b, minor semi-axis of ellipse
%       center, center of ellipse
%       theta, orientiation of ellipse
%       T, parametrisation [0,2*pi]
%
% Output:
%       z, function returning position at t
%       zp, function returning first derivative at t
%       zpp, function returning second derivative at t
%       
% Updated: - Anna June 3, 2020
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if nargin > 5
    if varargin{1} == 1
        dir = 1;
    else 
        dir = -1;
    end
else
    dir = 1;
end

A1 = (a+b)/2;
A2 = (a-b)/2;

if dir == 1
    z = (A1*exp(1i*T)+A2*exp(-1i*T))*exp(1i*theta)+center;
    zp = 1i*((A1*exp(1i*T)-A2*exp(-1i*T)))*exp(1i*theta);
    zpp = (-A1*exp(1i*T)-A2*exp(-1i*T))*exp(1i*theta);
else
    z = (A1*exp(-1i*T)+A2*exp(1i*T))*exp(1i*theta)+center;
    zp = 1i*((-A1*exp(-1i*T)+A2*exp(1i*T)))*exp(1i*theta);
    zpp = (-A1*exp(-1i*T)-A2*exp(1i*T))*exp(1i*theta);
end   
