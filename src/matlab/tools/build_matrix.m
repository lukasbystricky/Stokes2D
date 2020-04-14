function M = build_matrix(matvec, domain)

N = length(domain.z);
Q = eye(2*N);
M = zeros(2*N);

for i = 1:2*N
    M(:,i) = matvec(Q(:,i));
end