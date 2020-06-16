function m = matvec_double_layer_mobility(X, domain)
% MATVEC_DOUBLE_LAYER_MOBILITY: perform a matrix vector product using the 
% double-layer formulation to compute the mobility problem on an unbounded 
% domain, i.e. u_bg(x) = -q/2 + -D[q](x) - u_trans - omega(x - c)^\perp. 
% Here u_bg is the backgroudn flow, u_trans is the translational velocity,
% omega is the rotational velocity, and c is the center of the particles. 
% Uses FMM to compute the product quickly. Note we assume the particles are
% force- and torque-free.
%
% inputs:
% -X: vector containing [q1; q2; u_trans; omega]
% -domain: structure containing geometry information
%
% output:
% -m: vector containing [u1; u2; integral of q; integral of q(x-c)^perp]

m = zeros(length(X),1);
N = length(domain.z);
q = X(1:2*N);
n_particles = size(domain.wall_indices,1);

% extract translational and rotational velocities
start = 2*N + 1;
u_trans1 = X(start: start+n_particles-1);
u_trans2 = X(start+n_particles:start+2*n_particles - 1);
omega = X(start+2*n_particles: end);

x = real(domain.z);
y = imag(domain.z);
n1 = real(-1i*domain.zp)./abs(domain.zp);
n2 = imag(-1i*domain.zp)./abs(domain.zp);
q1 = q(1:N);
q2 = q(N+1:2*N);
quad_weights = domain.quad_weights;

% certain things are easier to do in complex variables
qc = q1 + 1i*q2;
qwazp = qc.*domain.wazp;
zp = domain.zp;
zpp = domain.zpp;

[u1, u2] = stokesDLPfmm(real(qwazp(:)),imag(qwazp(:)),x(:),y(:),n1(:),n2(:));
u = -u1 + -1i*u2;

%% correct using special quadrature
for i = 1:n_particles
    particle_indices = domain.wall_indices(i,1):domain.wall_indices(i,2);
    other_indices = setdiff(1 : length(qc), particle_indices);
    qc_particle = qc(particle_indices);
    u_target = u(other_indices);
    
    domainTmp = domain.particles{i};
    
    if isreal(u_target)
        u_target = u_target + 1e-300*1i;
    end
    
    if isreal(qc_particle)
        qc_particle = qc_particle + 1e-300*1i;
    end
    
    [u_new,~] = mex_SQ_dlp(domain.z(other_indices), domainTmp.z, domainTmp.zp,...
        domainTmp.quad_weights, domainTmp.panel_breaks, domainTmp.wazp, domainTmp.z32,...
        domainTmp.zp32, domainTmp.quad_weights32, domainTmp.wazp32,...
        qc_particle, u_target,domainTmp.mean_panel_length,...
        domainTmp.extra.gridSolidmat, ...
        domainTmp.extra.Nrows,domainTmp.extra.Ncols,domainTmp.extra.panels2wall,...
        domainTmp.reference_cell,false);
    
    u(other_indices) = u_new;
end

%% diagonal elements of the double-layer
u = u + quad_weights.*...
    (qc.*imag(zpp./zp) + ...
    conj(qc).*imag(zpp.*conj(zp))./conj(zp).^2)/(4*pi);

%% compute rotational velocity

rot1 = zeros(N,1);
rot2 = zeros(N,1);
u_trans1_vec = zeros(N,1);
u_trans2_vec = zeros(N,1);

for i = 1:n_particles
    particle_indices = domain.wall_indices(i,1):domain.wall_indices(i,2);
    
    x_wall = x(particle_indices);
    y_wall = y(particle_indices);
    
    c1 = real(domain.centers(i));
    c2 = imag(domain.centers(i));
    
    rot1(particle_indices) = -omega(i)*(y_wall - c2);
    rot2(particle_indices) = omega(i)*(x_wall - c1);
    
    u_trans1_vec(particle_indices) = u_trans1(i);
    u_trans2_vec(particle_indices) = u_trans2(i);
end

m(1:2*N) = [-q1/2 + real(u) - u_trans1_vec - rot1; -q2/2 + imag(u) - u_trans2_vec - rot2];


%% add density constraints to close system
for i = 1:n_particles
    particle_indices = domain.wall_indices(i,1):domain.wall_indices(i,2);
    
    qw_wall = qwazp(particle_indices);
    x_wall = x(particle_indices);
    y_wall = y(particle_indices);
    
    c1 = real(domain.centers(i));
    c2 = imag(domain.centers(i));
    
    m(2*N + i) = sum(real(qw_wall));
    m(2*N + n_particles + i) = sum(imag(qw_wall));
    m(2*N + 2*n_particles + i) = sum(-real(qw_wall).*(y_wall - c2) +...
                imag(qw_wall).*(x_wall - c1));
end

