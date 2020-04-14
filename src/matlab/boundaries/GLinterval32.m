function [t,w] = GLinterval32(a,b,NPanels,pedges)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [t,w] = GLinterval(a,b,NPanels)
%   Returns nodes and weights for composite 32 point Gauss-Legendre quadrature
%   using NPanels panels over (a,b)
%
% [t,w] = GLinterval(a,b,NPanels,pedges)
%   Returns nodes and weights for composite 32 point Gauss-Legendre quadrature using NPanels panels over 
%   (pedges(1),pedges(Npanels+1))  with the not necessarily equisized panels 
%   having panel edges in parameter as dictated in pedges.
%
% Input:
%       a,b, interval edges
%       Npanels, number of panels to divide interval into
%       pedges (optional), pre-set panel edges
%
% Output:
%       t, array of assembled 32 G-L nodes per panel
%       w, array of assembled 32 G-L weights per panel
%
% Updated: -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

[T,W] = GLinit32;

if nargin == 3

    t = zeros(32*NPanels,1);
    w = repmat(W*(b-a)/2/NPanels,NPanels,1);
    
    pedges = linspace(a,b,NPanels+1);
    
    ptr = 1;
    for j = 1:NPanels
        t(ptr:ptr+31) = (pedges(j+1)-pedges(j))*T/2+(pedges(j)+pedges(j+1))/2;
        ptr = ptr + 32;
    end

else
    t = zeros(32*NPanels,1);
    w = zeros(32*NPanels,1);
    
    ptr = 1;
    for j = 1:NPanels
        t(ptr:ptr+31) = (pedges(j+1)-pedges(j))*T/2+(pedges(j)+pedges(j+1))/2;
        w(ptr:ptr+31) = W*(pedges(j+1)-pedges(j))/2;
        ptr = ptr + 32;
    end
end


