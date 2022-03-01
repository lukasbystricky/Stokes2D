function m = matvec_double_layer_resistance_slip(X, problem)
% MATVEC_DOUBLE_LAYER_RESISTANCE_SLIP: perform a matrix vector product for
% the slip boundary condition  non-periodic problems. Uses the double-layer
% formulation with Power-Miranda completion, i.e. u(x) = -q/2 + -D[q](x) + 
% S[F](x) + R[T](x) where F is the net force on all inner walls, and T is 
% the net torque on all innerwalls. Uses FMM to compute this product
% quickly.
%
% inputs:
% -X: vector containing [q1; q2; F; T].
% -problem: problem structure.
%
% output:
% -m: vector containing [b2; b1; integral of q; integral of q(x-c)^perp],
%   where b1 is  (u \dot \nu) + alpha*(\nu \dot (\sigma \dot n)) and b1 is 
%   (u \dot \nu).

m = zeros(length(X),1);

domain = problem.domain;
x = real(domain.z);
y = imag(domain.z);
N = length(domain.z);

% extract net forces and torques
n_inner_walls = size(domain.wall_indices,1) - 1;
if n_inner_walls > 0
    start = 2*N + 1;
    F1 = X(start:start+n_inner_walls-1);
    F2 = X(start+n_inner_walls:start+2*n_inner_walls-1);
    T = X(start+2*n_inner_walls:end);
end

% density components
q = X(1:2*N);
q1 = q(1:N);
q2 = q(N+1:2*N);

% certain things are easier to do in complex variables
qc = q1 + 1i*q2;
qwazp = qc.*domain.wazp;

% normalized tangent
nu1 = real(domain.zp)./abs(domain.zp);
nu2 = imag(domain.zp)./abs(domain.zp);

% normalized normal
n1 = real(-1i*domain.zp)./abs(domain.zp);
n2 = imag(-1i*domain.zp)./abs(domain.zp);

% temporary solution structure
solution.problem = problem;
solution.q = [q1 q2];
solution.forces = F1 + 1i*F2;
solution.torques = T;
solution.local_indices = 1:length(solution.q);
solution.alpha = problem.alpha;

% velocity components (including jump)
[u1,u2] = evaluate_velocity_on_surface(solution, solution);

% normal velocity (u \dot n)
udotn = u1.*n1 + u2.*n2;
b1 = udotn;

% tangential velocity (u \dot \nu)
udotnu = u1.*nu1 + u2.*nu2;
b2 = udotnu;

if nnz(problem.alpha) > 0
    % stress components (including jump)
    [sxx,sxy,syx,syy] = evaluate_stress_on_surface(solution, solution, 'fluid');
    
    % traction t = (\sigma \dot n)
    t1 = n1.*sxx + n2.*syx;
    t2 = n1.*sxy + n2.*syy;

    % tangential component of traction (\nu \dot t)
    nudott = t1.*nu1 + t2.*nu2;
    
    % check alpha
    numalpha = -mean(real(udotnu))/mean(nudott);

    % add second term
    b2 = b2 + problem.alpha.*nudott;
    
    fprintf('udotn=%f,    udotnu=%f,    nudott=%f,    alpha=%f, nudott*alpha=%f\n', ...
        mean(udotn(1:end/2)), mean(udotnu(1:end/2)), mean(nudott(1:end/2)), -mean(real(udotnu(1:end/2)))/mean(nudott(1:end/2)), mean(nudott(1:end/2))*numalpha);
    fprintf('      %f,           %f,           %f,          %f,              %f\n', ...
        mean(udotn(end/2+1:end)), mean(udotnu(end/2+1:end)), mean(nudott(end/2+1:end)), -mean(real(udotnu(end/2+1:end)))/mean(nudott(end/2+1:end)), mean(nudott(end/2+1:end))*numalpha);
else
    fprintf('udotn=%f,    udotnu=%f\n', mean(udotn(1:end/2)), mean(udotnu(1:end/2)));
    fprintf('      %f,           %f\n', mean(udotn(end/2+1:end)), mean(udotnu(end/2+1:end)));
end

% add Robin boundary condition
m(1:2*N) = [b2; b1];

% add density constraints to close system
for i = 1:n_inner_walls
    qw_wall = qwazp(domain.wall_indices(i+1,1):domain.wall_indices(i+1,2));
    x_wall = x(domain.wall_indices(i+1,1):domain.wall_indices(i+1,2));
    y_wall = y(domain.wall_indices(i+1,1):domain.wall_indices(i+1,2));

    c1 = real(domain.centers(i+1));
    c2 = imag(domain.centers(i+1));

    m(2*N+i) = sum(real(qw_wall)) - F1(i);
    m(2*N+n_inner_walls+i) = sum(imag(qw_wall)) - F2(i);
    m(2*N+2*n_inner_walls+i) = sum(-real(qw_wall).*(y_wall - c2) + ...
        imag(qw_wall).*(x_wall - c1)) - T(i);
end
