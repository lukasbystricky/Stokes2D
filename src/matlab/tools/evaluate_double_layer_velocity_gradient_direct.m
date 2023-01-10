function dV = evaluate_double_layer_velocity_gradient_direct(xtar, ytar, z, btar, zp, q, w)
ztar = xtar + 1i*ytar;
n = -1i*zp./abs(zp);

% inefficient, but will do for now
Ih_q = evaluate_hypersingular_integral(xtar, ytar, z, zp, q, w);
Ih_rqn = evaluate_hypersingular_integral(xtar, ytar, z, zp, real(q.*conj(n))./n, w);
Is_q = evaluate_supersingular_integral(xtar, ytar, z, zp, q, w);
Is_tq = evaluate_supersingular_integral(xtar, ytar, z, zp, conj(z).*q, w);

dV = -1i/(2*pi) * (btar.*real(Ih_q) - conj(btar.*(Ih_rqn+conj(ztar).*Is_q-Is_tq)));
end

function  Ih = evaluate_hypersingular_integral(xtar, ytar, z, zp, q, w)
Ih = zeros(size(xtar));
for k = 1:length(xtar)
    rho = (xtar(k) + 1i*ytar(k)) - z;
    Ih(k) = sum(q.*zp.*w./rho.^2);
end
end

function  Is = evaluate_supersingular_integral(xtar, ytar, z, zp, q, w)
Is = zeros(size(xtar));
for k = 1:length(xtar)
    rho = (xtar(k) + 1i*ytar(k)) - z;
    Is(k) = sum(q.*zp.*w./rho.^3);
end
end