function M = build_matrix(matvec, domain)

N = 2*length(domain.z) + 6;
Q = eye(N);
M = zeros(N);

for i = 1:N
    M(:,i) = matvec(Q(:,i));
end