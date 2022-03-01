function sigma = evaluate_double_layer_stress_direct(xtar, ytar, z, btar, zp, q, w)
ztar = xtar + 1i*ytar;
n = -1i*zp./abs(zp);

Ih1 = evaluate_hypersingular_integral(xtar, ytar, z, zp, q, w);
Ih2 = evaluate_hypersingular_integral(xtar, ytar, z, zp, real(q.*conj(n))./n, w);
Is1 = evaluate_supersingular_integral(xtar, ytar, z, zp, conj(z).*q, w);
Is2 = evaluate_supersingular_integral(xtar, ytar, z, zp, q, w);

sigma = -(btar.*imag(Ih1)/2 - 1i*conj(btar.*(Ih2+Is1-conj(ztar).*Is2)))/pi;

end

function  Ih = evaluate_hypersingular_integral(xtar, ytar, z, zp, q, w)
Ih = zeros(size(xtar));
indices = 1:length(xtar);
for k = 1:length(xtar)
    indices_tmp = indices;
    indices_tmp(indices==k) = [];
    rho = (xtar(k) + 1i*ytar(k)) - z(indices_tmp);
    Ih(k) = sum(q(indices_tmp).*zp(indices_tmp).*w(indices_tmp)./rho.^2);
end
end

function  Is = evaluate_supersingular_integral(xtar, ytar, z, zp, q, w)
Is = zeros(size(xtar));
indices = 1:length(xtar);
for k = 1:length(xtar)
    indices_tmp = indices;
    indices_tmp(indices==k) = [];
    rho = (xtar(k) + 1i*ytar(k)) - z(indices_tmp);
    Is(k) = sum(q(indices_tmp).*zp(indices_tmp).*w(indices_tmp)./(-rho).^3);
end
end