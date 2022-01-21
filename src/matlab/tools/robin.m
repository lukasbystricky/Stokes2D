function rhs = robin(X, problem)
% ROBIN: sets up the right hand side when using Robin boundary conditions
% inputs:
%
% -X: vector containing [q1; q2; <u1>; <u2>]
% -problem: problem strucutre
%
% output:
% -rhs: vector containing [b1; b2; 0; 0], where b1 is (u \dot n)
% and b2 is (u \dot \nu) + alpha*(\nu \dot (\sigma \dot n)).

domain = problem.domain;
N = length(domain.z);

% density components
q = X(1:2*N);
q1 = q(1:end/2);
q2 = q(end/2+1:end);

% normalized tangent
nu1 = real(domain.zp)./abs(domain.zp);
nu2 = imag(domain.zp)./abs(domain.zp);    

% normalized normal
n1 = real(-1i*domain.zp)./abs(domain.zp);
n2 = imag(-1i*domain.zp)./abs(domain.zp);

% temporary solution structure
solution.problem = problem;
solution.q = [q1 q2];
solution.u_avg = [X(end-1); X(end)];
solution.local_indices = 1:length(solution.q);

[u1,u2] = evaluate_velocity_on_surface(solution, solution);

% normal velocity (u \dot n)
b1 = u1.*n1 + u2.*n2;

% tangential velocity (u \dot \nu)
b2 = u1.*nu1 + u2.*nu2;

if problem.alpha > 0
    [sxx,sxy,syx,syy] = evaluate_stress_on_surface(solution, solution, 'fluid');

    % traction t = (\sigma \dot n)
    t1 = n1.*sxx + n2.*sxy;
    t2 = n1.*syx + n2.*syy;

    % tangential component of traction times alpha (alpha * (\nu \dot t))
    tnu = t1.*nu1 + t2.*nu2;
    
    % check alpha
    numalpha = -mean(real(b2))/mean(tnu)
    
    % add second term
    b2 = b2 + problem.alpha*tnu;
end

% fix dimension by appending two zeros due to pressure constraint
rhs = [b1; b2; 0; 0];
