function [duS, duR] = completion_contribution_velocity_gradient(zsrc, ztar, btar, forces, torques)

duS = zeros(length(ztar),1);
duR = zeros(length(ztar),1);

for i = 1:length(zsrc)
    r = ztar - zsrc(i);
    
    duS = duS - ((r*conj(forces(i)).*conj(btar))./(2*conj(r.^2)) + ...
        (forces(i)*conj(btar)-btar*conj(forces(i)))./(2*conj(r)) + ...
        (btar*forces(i))./(2*r));

    duR = duR - 1i*2*torques(i)*conj(btar./r.^2)/(4*pi);
end