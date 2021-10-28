function omega = evaluate_double_layer_vorticity_direct(xtar, ytar, z, zp, q, w)
Ih_f = evaluate_hypersingular_integral(xtar, ytar, z, zp, q, w);
omega = -real(Ih_f)/pi;
end

function  Ih = evaluate_hypersingular_integral(xtar, ytar, z, zp, q, w)
Ih = zeros(size(xtar));
for k = 1:length(xtar)
    rho = (xtar(k) + 1i*ytar(k)) - z;
    Ih(k) = sum(q.*zp.*w./rho.^2);
end
end