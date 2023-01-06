function P = evaluate_double_layer_pressure_direct(xtar, ytar, z, zp, q, w)

P = zeros(size(xtar));
for k = 1:length(xtar)

    rho = (xtar(k) + 1i*ytar(k)) - z;
    P(k) = imag(sum(q.*zp.*w./(rho.^2)))/pi;
end


