function [t,w,pedges] = GLinterval(a,b,NPanels,pedges)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [t,w,pedges] = GLinterval(a,b,NPanels)
%   Returns nodes, weights and panel edges for composite 16 point 
%   Gauss-Legendre quadrature using NPanels equisized panels over (a,b)
%
% [t,w,pedges] = GLinterval(a,b,NPanels,pedges)
%   Returns nodes, weights and panel edges for composite 16 point 
%   Gauss-Legendre quadrature using NPanels panels over 
%   (pedges(1),pedges(Npanels+1))  with the not necessarily equisized panels 
%   having panel edges in parameter as dictated in pedges.
%
% Input:
%       a,b, interval edges
%       Npanels, number of panels to divide interval into
%       pedges (optional), pre-set edges of panels (not necessarily
%       equisized)
% 
% Output:
%       t, array with 16 G-L nodes per panel
%       w, array with 16 G-L weights per panel
%       pedges, panel edges
%
% Updated: -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

[T,W] = GLinit16;

if nargin == 3

    t = zeros(16*NPanels,1);
    w = repmat(W*(b-a)/2/NPanels,NPanels,1);
    
    pedges = linspace(a,b,NPanels+1);
    
    ptr = 1;
    for j = 1:NPanels
        t(ptr:ptr+15) = (pedges(j+1)-pedges(j))*T/2+(pedges(j)+pedges(j+1))/2;
        ptr = ptr + 16;
    end

else
    t = zeros(16*NPanels,1);
    w = zeros(16*NPanels,1);
    
    ptr = 1;
    for j = 1:NPanels
        t(ptr:ptr+15) = (pedges(j+1)-pedges(j))*T/2+(pedges(j)+pedges(j+1))/2;
        w(ptr:ptr+15) = W*(pedges(j+1)-pedges(j))/2;
        ptr = ptr + 16;
    end
end

