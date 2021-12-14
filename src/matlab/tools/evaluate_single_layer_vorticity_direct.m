function omega = evaluate_single_layer_vorticity_direct(xtar, ytar, z, zp, q, w)
n = -1i*zp./abs(zp);
Ic = evaluate_cauchy_integral(xtar, ytar, z, zp, q./n, w);
omega = -real(Ic)/(2*pi);
end

function  Ic = evaluate_cauchy_integral(xtar, ytar, z, zp, q, w)
Ic = zeros(size(xtar));
for k = 1:length(xtar)
    rho = (xtar(k) + 1i*ytar(k)) - z;
    Ic(k) = sum(q.*zp.*w./-rho);
end
end