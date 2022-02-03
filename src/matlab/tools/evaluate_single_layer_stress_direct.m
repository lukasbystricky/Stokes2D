function sigma = evaluate_single_layer_stress_direct(xtar, ytar, z, btar, zp, q, w)
ztar = xtar + 1i*ytar;
n = -1i*zp./abs(zp);

Ic1 = evaluate_cauchy_integral(xtar, ytar, z, zp, q./n, w);
Ic2 = evaluate_cauchy_integral(xtar, ytar, z, zp, conj(q)./n, w);
Ih1 = evaluate_hypersingular_integral(xtar, ytar, z, zp, conj(z.*q)./n, w);
Ih2 = evaluate_hypersingular_integral(xtar, ytar, z, zp, conj(q)./n, w);

sigma = (2*btar.*imag(Ic1) + 1i*conj(btar.*(Ic2+Ih1-conj(ztar).*Ih2)))/(4*pi);

end

function  Ic = evaluate_cauchy_integral(xtar, ytar, z, zp, q, w)
Ic = zeros(size(xtar));
for k = 1:length(xtar)
    rho = -((xtar(k) + 1i*ytar(k)) - z);
    Ic(k) = sum(q.*zp.*w./rho);
end
end

function  Ih = evaluate_hypersingular_integral(xtar, ytar, z, zp, q, w)
Ih = zeros(size(xtar));
for k = 1:length(xtar)
    rho = (xtar(k) + 1i*ytar(k)) - z;
    Ih(k) = sum(q.*zp.*w./rho.^2);
end
end