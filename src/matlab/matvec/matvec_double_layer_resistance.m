function m = matvec_double_layer_resistance(X, domain, fmm)
% MATVEC_DOUBLE_LAYER: perform a matrix vector product using the 
% double-layer formulation with Power-Miranda completion, i.e. 
% u(x) = -q/2 + -D[q](x) + S[F](x) + R[T](x), F is the net force on all 
% inner walls, and T is the net torque on all inner walls. Uses FMM to
% compute this product quickly.
%
% inputs:
% -X: vector containing [q1; q2; F; T]
% -domain: structure containing geometry information
%
% output:
% -m: vector containing [u1; u2; integral of q; integral of q(x-c)^perp]

m = zeros(length(X),1);
N = length(domain.z);
q = X(1:2*N);

% extract net forces and torques
n_inner_walls = domain.n_inner_walls;

if n_inner_walls > 0
    start = 2*N + 1;
    F1 = X(start: start+n_inner_walls-1);
    F2 = X(start+n_inner_walls:start+2*n_inner_walls - 1);
    T = X(start+2*n_inner_walls: end);
end

x = real(domain.z);
y = imag(domain.z);
n1 = real(-1i*domain.zp)./abs(domain.zp);
n2 = imag(-1i*domain.zp)./abs(domain.zp);
q1 = q(1:N);
q2 = q(N+1:2*N);

% certain things are easier to do in complex variables
qc = q1 + 1i*q2;
qwazp = qc.*domain.wazp;

outer_wall_indices = domain.outer_wall_indices;
inner_wall_indices = domain.inner_wall_indices;

if fmm
    [u1, u2] = stokesDLPfmm(real(qwazp(:)),imag(qwazp(:)),x(:),y(:),n1(:),n2(:));
else
    u1 = zeros(numel(q1),1);
    u2 = zeros(numel(q1),1);
    for k = 1:numel(q1)
        
        % skip self interaction term
        ind = [(1:k-1) (k+1:numel(q1))];
        
        rx = x(k) - x(ind);
        ry = y(k) - y(ind);
        rho4 = (rx.^2 + ry.^2).^2;
        
        rdotq = rx.*real(qwazp(ind)) + ry.*imag(qwazp(ind));
        rdotn = rx.*n1(ind) + ry.*n2(ind);
        
        u1(k) = 4*sum(rdotn.*rdotq./rho4.*rx);
        u2(k) = 4*sum(rdotn.*rdotq./rho4.*ry);
    end
    
    u1 = u1/4/pi;
    u2 = u2/4/pi;
end

u = -u1 + -1i*u2;

%% diagonal elements of the double-layer
% NB: Not 100% sure why we have to split into the two cases here and not in
% periodic case
quad_weights_in = domain.quad_weights(inner_wall_indices);
zp_in = domain.zp(inner_wall_indices);
zpp_in = domain.zpp(inner_wall_indices);
qc_in = qc(inner_wall_indices);
quad_weights_out = domain.quad_weights(outer_wall_indices);
zp_out = domain.zp(outer_wall_indices);
zpp_out = domain.zpp(outer_wall_indices);
qc_out = qc(outer_wall_indices);

u(inner_wall_indices) = u(inner_wall_indices) + quad_weights_in.*...
    (qc_in.*imag(zpp_in./zp_in) + ...
    conj(qc_in).*imag(zpp_in.*conj(zp_in))./conj(zp_in).^2)/(4*pi);
u(outer_wall_indices) = u(outer_wall_indices) - quad_weights_out.*...
    (qc_out.*imag(zpp_out./zp_out) + ...
    conj(qc_out).*imag(zpp_out.*conj(zp_out))./conj(zp_out).^2)/(4*pi);

%% remove container nullspace
n_outer1 = n1(outer_wall_indices);
n_outer2 = n2(outer_wall_indices);
qw_outer1 = real(qwazp(outer_wall_indices));
qw_outer2 = imag(qwazp(outer_wall_indices));
N0 = (n_outer1 + 1i*n_outer2).*(n_outer1.*qw_outer1 + n_outer2.*qw_outer2);

% This seems to be incorrect, which is why it is commented out. DK
%u(outer_wall_indices) = u(outer_wall_indices) + N0;

%% add rotlet and Stokelet contribution and put everything back to real 
% variables
if n_inner_walls > 0
    [uS, uR] = completion_contribution(domain.centers(2:end), x + 1i*y,...
        F1 + 1i*F2, T);
else
    uS = 0;
    uR = 0;
end

m(1:2*N) = [-q1/2 + real(u + uS + uR); -q2/2 + imag(u + uS + uR)];

%% add density constraints to close system
for i = 1:n_inner_walls
    
    qw_wall = qwazp(domain.wall_indices(i+1,1):domain.wall_indices(i+1,2));
    x_wall = x(domain.wall_indices(i+1,1):domain.wall_indices(i+1,2));
    y_wall = y(domain.wall_indices(i+1,1):domain.wall_indices(i+1,2));
    
    c1 = real(domain.centers(i+1));
    c2 = imag(domain.centers(i+1));
    
    m(2*N + i) = sum(real(qw_wall)) - F1(i);
    m(2*N + n_inner_walls + i) = sum(imag(qw_wall)) - F2(i);
    m(2*N + 2*n_inner_walls + i) = sum(-real(qw_wall).*(y_wall - c2) +...
                imag(qw_wall).*(x_wall - c1)) - T(i);
end

