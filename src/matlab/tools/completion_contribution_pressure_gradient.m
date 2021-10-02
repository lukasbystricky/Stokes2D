function dPs = completion_contribution_pressure_gradient(zsrc, ztar, forces)
dPs = zeros(length(ztar),1);
for i = 1:length(zsrc)
    rho = ztar - zsrc(i);
    dPs = dPs - 1/(2*pi)*conj(forces(i)./rho.^2);
end