function [uS, uR] = completion_contribution(zsrc, ztar, forces, torques)

uS = zeros(length(ztar),1);
uR = zeros(length(ztar),1);

for i = 1:length(zsrc)
    r = ztar - zsrc(i);
    
    uS = uS + -log(abs(r))*forces(i) + r.*(real(forces(i))*real(r) +...
                        imag(forces(i))*imag(r))./abs(r).^2;
                    
    uR = uR + (imag(r) - 1i*real(r))*torques(i)./abs(r).^2;
end

