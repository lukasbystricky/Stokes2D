function omegaS = completion_contribution_vorticity(zsrc, ztar, forces)
omegaS = zeros(length(ztar),1);

for i = 1:length(zsrc)
    r = ztar - zsrc(i);
    
    rotlet_kernel = 1i*r./abs(r).^2;
    omegaS = omegaS + 2*rotlet_kernel*forces(i)/(4*pi);
end