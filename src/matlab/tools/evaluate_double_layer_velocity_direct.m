function u = evaluate_double_layer_velocity_direct(xtar, ytar, z, zp, q, w)
ztar = xtar + 1i*ytar;
n = -1i*zp./abs(zp);

Ic1 = evaluate_cauchy_integral(xtar, ytar, z, zp, q, w);
Ih1 = evaluate_hypersingular_integral(xtar, ytar, z, zp, q, w);
Ih2 = evaluate_hypersingular_integral(xtar, ytar, z, zp, conj(z).*q, w);
Ic2 = evaluate_cauchy_integral(xtar, ytar, z, zp, real(q.*conj(n)), w);

u = 1i*(Ic1 + conj(conj(ztar).*Ih1-Ih2-2*Ic2))/(4*pi);

end

function  Is = evaluate_cauchy_integral(xtar, ytar, z, zp, q, w)
Is = zeros(size(xtar));
for k = 1:length(xtar)
    rho = (xtar(k) + 1i*ytar(k)) - z;
    Is(k) = sum(q.*zp.*w./(-rho));
end
end

function  Ih = evaluate_hypersingular_integral(xtar, ytar, z, zp, q, w)
Ih = zeros(size(xtar));
for k = 1:length(xtar)
    rho = (xtar(k) + 1i*ytar(k)) - z;
    Ih(k) = sum(q.*zp.*w./rho.^2);
end
end
