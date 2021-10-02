function dP = evaluate_double_layer_pressure_gradient_direct(xtar, ytar, z, zp, q, w)
dP = -(2*1i)/pi * conj(evaluate_supersingular_integral(xtar,ytar,z,zp,q,w));
end

function  Is = evaluate_supersingular_integral(xtar, ytar, z, zp, q, w)
Is = zeros(size(xtar));
for k = 1:length(xtar)
    r = z - (xtar(k) + 1i*ytar(k));
    Is(k) = sum(q.*zp.*w./r.^3);
end
end