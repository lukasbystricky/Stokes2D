function [duS, duR] = completion_contribution_gradient(zsrc, ztar, btar, forces, torques)

duS = zeros(length(ztar),1);
duR = zeros(length(ztar),1);

for i = 1:length(zsrc)
    r = ztar - zsrc(i);
    
    %uS = uS + (-log(abs(r))*forces(i) + real(r*conj(forces(i)))./conj(r)) /4/pi;
    duS = duS - real((r*conj(forces(i)).*conj(btar))./(2*conj(r.^2)) + ...
        (forces(i)*conj(btar)-btar*conj(forces(i)))./(2*conj(r)) + ...
        (btar*forces(i))./(2*r));
    
    %uR = uR + (-imag(r) + 1i*real(r))*torques(i)./abs(r).^2/2/pi;
    %duR = duR - 1/(2*pi) * real((torques(i)*btar(i))./r.^2 + ...
    %    (conj(torques(i))*conj(btar))./conj(r.^2));
    duR = duR - torques(i)./conj(r.^2);
end