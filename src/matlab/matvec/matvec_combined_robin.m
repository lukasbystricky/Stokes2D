function m = matvec_combined_robin(X, problem)
% MATVEC_COMBINED: perform a matrix vector product using the combined-layer
% formulation, i.e. u(x) = q/2 + eta*S[q](x) + D[q](x) + <u>, where <u> is 
% the average velocity over a reference cell. Uses periodic spectral Ewald 
% to compute the product quickly.
%
% inputs:
% -X: vector containing [q1; q2; <u1>; <u2>]
% -domain: structure containing geometry information
% -eta: scaling parameter in front of single-layer potential
%
% output:
% -m: vector containing [u1; u2; dp1; dp2], where u1 and u2 are the 
%   velocity component on the boundary, and dp1 and dp2 are the average 
%   pressure gradients in the x and y direction

domain = problem.domain;
eta = problem.eta;

m = zeros(length(X),1);
N = length(domain.z);

q = X(1:2*N);
u_avg = X(end-1:end);
V = domain.Lx * domain.Ly;

x = real(domain.z);
y = imag(domain.z);
n1 = real(-1i*domain.zp)./abs(domain.zp);
n2 = imag(-1i*domain.zp)./abs(domain.zp);
q1 = q(1:end/2);
q2 = q(end/2+1:end);

% certain things are easier to do in complex variables
qc = q1 + 1i*q2;
qwazp = qc.*domain.wazp;

[u1, u2] = StokesDLP_ewald_2p(x, y, x, y, n1, n2, real(qwazp), imag(qwazp),...
    domain.Lx, domain.Ly);

if ~isinf(eta)
    T = u1 + 1i*u2;
    
    % diagonal elements of the double-layer
    T = T + domain.quad_weights.*(qc.*imag(domain.zpp./domain.zp) + ...
        conj(qc).*imag(domain.zpp.*conj(domain.zp))./conj(domain.zp).^2)/(4*pi);
end

[u1, u2] = StokesSLP_ewald_2p(x, y, x, y, real(qwazp), imag(qwazp),...
    domain.Lx, domain.Ly);
G = u1 + 1i*u2;

% apply log-quadrature corrections for the Stokeslet
wall_indices = domain.wall_indices;
for i = 1:size(wall_indices,1)
    
    panels_per_wall = (wall_indices(i,2)-wall_indices(i,1)+1)/16;
    thiswall = wall_indices(i,1):wall_indices(i,2);
    G(thiswall) = G(thiswall) - qwazp(thiswall).*...
        log(pi/panels_per_wall*abs(domain.zp(thiswall)))/(4*pi);
    
    for j = 1:panels_per_wall
        
        source_ind = wall_indices(i,1)-1+((j-1)*16+(1:16));
        target_ind = wall_indices(i,1)+...
                mod((j-1)*16+(-3:20)+panels_per_wall*16-1, panels_per_wall*16);
        
        correction = -domain.Lmod*(qwazp(source_ind)) / (4*pi);
        
        G(target_ind) = G(target_ind) + correction;
    end
end

G = G + (qc + domain.zp./conj(domain.zp).*conj(qc)).*domain.wazp/(8*pi);

if ~isinf(eta)
    u = -qc/2 + eta*G + T;
else
    u = G;
end

% put everything back to real variables
m(1:2*N) = [real(u) + u_avg(1); imag(u) + u_avg(2)];

% add pressure constraint
if ~isinf(eta)
    m(end - 1) = eta*real(sum(qwazp)) / V;
    m(end) = eta*imag(sum(qwazp)) / V;
else
    m(end - 1) = real(sum(qwazp)) / V;
    m(end) = imag(sum(qwazp)) / V;
end

% Robin boundary condition
bc = robin(X, problem);
m = m + bc;
