function sigma = evaluate_double_layer_stress_direct(xtar, ytar, z, btar, zp, q, w)
ztar = xtar + 1i*ytar;
n = -1i*zp./abs(zp);

Ih1 = evaluate_hypersingular_integral(xtar, ytar, z, zp, q, w);
Ih2 = evaluate_hypersingular_integral(xtar, ytar, z, zp, real(q.*conj(n))./n, w);
Is1 = evaluate_supersingular_integral(xtar, ytar, z, zp, conj(z).*q, w);
Is2 = evaluate_supersingular_integral(xtar, ytar, z, zp, q, w);

sigma = -(btar.*imag(Ih1) + 1i*conj(btar.*(Ih2+Is1-conj(ztar).*Is2)))/pi;

% % test steps in derivation
% sigma2 = zeros(size(xtar));
% for j = 1:length(xtar)
%     rho = (xtar(j) + 1i*ytar(j)) - z;
%     r2 = rho.*conj(rho);
%     r4 = r2.^2;
%     r6 = r2.^3;
%     
%     qdotn = real(q.*conj(n));
%     rdotb = real(rho.*conj(btar(j)));
%     rdotn = real(rho.*conj(n));
%     rdotq = real(rho.*conj(q));
%     bdotn = real(btar(j).*conj(n));
%     qdotb = real(q.*conj(btar(j)));
% 
%     % real
%     sigma2(j) = sum(((btar(j)*qdotn)./r2 + (q.*rdotb.*rdotn+n.*rdotq.*rdotb)./r4 + ...
%         (rho.*rdotn.*qdotb+rho.*rdotq.*bdotn)./r4 - ...
%         8*(rho.*rdotb.*rdotq.*rdotn)./r6).*w.*zp./(1i*n));
% 
%     % replace 1st term
%     sigma2(j) = sum((btar(j)*(q.*conj(n)+conj(q).*n)./(2*rho.*conj(rho)) + (q.*rdotb.*rdotn+n.*rdotq.*rdotb)./r4 + ...
%         (rho.*rdotn.*qdotb+rho.*rdotq.*bdotn)./r4 - ...
%         8*(rho.*rdotb.*rdotq.*rdotn)./r6).*w.*zp./(1i*n));
%     
%     % replace 2nd term
%     sigma2(j) = sum((btar(j)*(q.*conj(n)+conj(q).*n)./(2*rho.*conj(rho)) + ...
%         0.5*((conj(btar(j))*q.*conj(n)+conj(btar(j))*conj(q).*n+btar(j)*conj(q.*n))./conj(rho).^2 + ...
%         (2*conj(btar(j))*q.*n+btar(j)*q.*conj(n)+btar(j)*conj(q).*n)./(rho.*conj(rho)) + ...
%         (btar(j)*q.*n)./rho.^2) - ...
%         8*(rho.*rdotb.*rdotq.*rdotn)./r6).*w.*zp./(1i*n));
%     
%     % replace 3rd term
%     sigma2(j) = sum((btar(j)*(q.*conj(n)+conj(q).*n)./(2*rho.*conj(rho)) + ...
%         0.5*((conj(btar(j))*q.*conj(n)+conj(btar(j))*conj(q).*n+btar(j)*conj(q.*n))./conj(rho).^2 + ...
%         (2*conj(btar(j))*q.*n+btar(j)*q.*conj(n)+btar(j)*conj(q).*n)./(rho.*conj(rho)) + ...
%         (btar(j)*q.*n)./rho.^2) - ...
%         (rho.*conj(btar(j)*q.*n))./conj(rho).^3 - ...
%         (conj(btar(j)*q).*n+conj(btar(j))*q.*conj(n)+btar(j)*conj(q.*n))./conj(rho).^2 - ...
%         (conj(btar(j))*q.*n+btar(j)*conj(q).*n+btar(j)*q.*conj(n))./(rho.*conj(rho)) - ...
%         (btar(j)*q.*n)./rho.^2).*w.*zp./(1i*n));
%     
%     % simplified
%     sigma2(j) = -0.5*sum(((conj(btar(j))*q.*conj(n)+conj(btar(j)*q).*n+btar(j)*conj(q.*n))./conj(rho).^2 + ... 
%         (btar(j)*q.*n)./rho.^2 + ...
%         2*(rho.*conj(btar(j)*q.*n))./conj(rho).^3).*w.*zp./(1i*n));
% end
% sigma2 = sigma2/pi;
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
    Is(k) = sum(q.*zp.*w./(-rho).^3);
end
end