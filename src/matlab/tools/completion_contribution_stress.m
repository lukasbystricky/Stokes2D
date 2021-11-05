function [duS, duR] = completion_contribution_stress(zsrc, ztar, btar, forces, torques)

duS = zeros(length(ztar),1);
duR = zeros(length(ztar),1);

for i = 1:length(zsrc)
    r = ztar - zsrc(i);
    
    duS = duS - (conj(btar*forces(i))./conj(r).^3 + ... 
        (conj(btar)*forces(i)+btar*conj(forces(i)))./(r.*conj(r).^2) + ...
        (btar*forces(i))./(r.^2.*conj(r)))/(4*pi);
    
    %duR = duR - (1i*2*torques(i)) * (conj(btar./r.^2) + btar./abs(r))/(2*pi);
    duR = duR - 1i*conj(btar./r.^2)*torques(i)/pi;
end