function omega = evaluate_single_layer_vorticity_direct(xtar, ytar, z, zp, q, w)
n = -1i*zp./abs(zp);

Ic = evaluate_cauchy_integral(xtar, ytar, z, zp, q./n, w);

omega = -real(Ic)/(2*pi);

% omega1 = zeros(size(xtar));
% omega2 = zeros(size(xtar));
% qperp = -1i*q;
% for j = 1:length(xtar)
%     rho = (xtar(j) + 1i*ytar(j)) - z;
%     rhoperp = -1i*rho;
%     r2 = rho.*conj(rho);
%     
%     qdotrperp = real(q.*conj(rhoperp));
%     rdotqperp = real(rho.*conj(qperp));
%     
%     omega1(j) = sum(((qdotrperp-rdotqperp)./r2).*w.*zp./(1i*n));
% end
% omega1 = omega1/(4*pi);
% norm(omega-omega1)
% omega;
end

function  Ic = evaluate_cauchy_integral(xtar, ytar, z, zp, q, w)
Ic = zeros(size(xtar));
for k = 1:length(xtar)
    rho = (xtar(k) + 1i*ytar(k)) - z;
    Ic(k) = sum(q.*zp.*w./-rho);
end
end