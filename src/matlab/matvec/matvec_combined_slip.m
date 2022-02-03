function m = matvec_combined_slip(X, problem)
% MATVEC_COMBINED_SLIP: perform a matrix vector product using the combined
% -layer formulation, i.e. u(x) = q/2 + eta*S[q](x) + D[q](x) + <u>, where
% <u> is the average velocity over a reference cell. Uses periodic spectral
% Ewald to compute the product quickly.
%
% inputs:
% -X: vector containing [q1; q2; <u1>; <u2>]
% -domain: structure containing geometry information
% -eta: scaling parameter in front of single-layer potential
%
% output:
% -m: vector containing [b1; b2; dp1; dp2], where u1 and u2 are the 
%   velocity component on the boundary, and dp1 and dp2 are the average 
%   pressure gradients in the x and y direction

% ROBIN: sets up the Robin boundary conditions
%
% inputs:
% -X: vector containing [q1; q2; <u1>; <u2>]
% -problem: problem strucutre
%
% output:
% -b: vector containing [b2; b1], where b1 is (u \dot \nu) + alpha*(\nu
%   \dot (\sigma \dot n)) and b1 is (u \dot \nu).

m = zeros(length(X),1);

domain = problem.domain;
eta = problem.eta;
N = length(domain.z);
V = domain.Lx * domain.Ly;

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
solution.u_avg = X(2*N+1:2*N+2);
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

%if sum(problem.alpha) ~= 0
if nnz(problem.alpha) > 0
    [sxx,sxy,syx,syy] = evaluate_stress_on_surface(solution, solution, 'fluid');
    
    % traction t = (\sigma \dot n)
    t1 = n1.*sxx + n2.*syx;
    t2 = n1.*sxy + n2.*syy;

    % tangential component of traction (\nu \dot t)
    nudott = t1.*nu1 + t2.*nu2;
    
%     [ux,uy,vx,vy] = evaluate_velocity_gradient_on_surface(solution, solution, 'fluid');
%     nudott = (vx+uy).*n2.*nu1;
    
    % check alpha
    numalpha = -mean(real(udotnu))/mean(nudott);

    % add second term
    b2 = b2 + problem.alpha.*nudott;
    
    fprintf('udotn=%f,    udotnu=%f,    nudott=%f,    alpha=%f, nudott*alpha=%f\n', ...
        mean(udotn(1:end/2)), mean(udotnu(1:end/2)), mean(nudott(1:end/2)), numalpha, mean(nudott(1:end/2))*numalpha);
    fprintf('      %f,           %f,           %f,          %f,              %f\n', ...
        mean(udotn(end/2+1:end)), mean(udotnu(end/2+1:end)), mean(nudott(end/2+1:end)), numalpha, mean(nudott(end/2+1:end))*numalpha);
else
    fprintf('udotn=%f,    udotnu=%f\n', mean(udotn(1:end/2)), mean(udotnu(1:end/2)));
    fprintf('      %f,           %f\n', mean(udotn(end/2+1:end)), mean(udotnu(end/2+1:end)));
end

% add Robin boundary condition
%m(1:2*N) = [b1; b2];   % worse
m(1:2*N) = [b2; b1];    % better

% add pressure constraint
if ~isinf(eta)
    m(2*N+1) = eta*real(sum(qwazp)) / V;
    m(2*N+2) = eta*imag(sum(qwazp)) / V;
else
    m(2*N+1) = real(sum(qwazp)) / V;
    m(2*N+2) = imag(sum(qwazp)) / V;
end
