function dV = evaluate_single_layer_velocity_gradient_direct(xtar, ytar, z, btar, zp, q, w)
ztar = xtar + 1i*ytar;
n = -1i*zp./abs(zp);

Ih_qn = evaluate_hypersingular_integral(xtar, ytar, z, zp, q./n, w);
Ih_tqn = evaluate_hypersingular_integral(xtar, ytar, z, zp, conj(z).*q./n, w);
Ic_cqn = evaluate_cauchy_integral(xtar, ytar, z, zp, conj(q)./n, w);
Ic_qn = evaluate_cauchy_integral(xtar, ytar, z, zp, q./n, w);

dV = -1i*conj(btar)/2 .* (ztar.*conj(Ih_qn) - conj(Ih_tqn) + conj(Ic_cqn)) + ...
    1i*btar.*real(Ic_qn);
dV = 1/(4*pi)*dV;
end

function  Ic = evaluate_cauchy_integral(xtar, ytar, z, zp, q, w)
Ic = zeros(size(xtar));
for k = 1:length(xtar)
    rho = (xtar(k) + 1i*ytar(k)) - z;
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