function [L] = create_Lmod(varargin)
% Creates the full matrix Lmod with corrections of the quadrature weights
% for computing the log-kernel. 

% Set target points (optional)
% ztar = [];

global high_precision;

high_precision = 1;

if nargin < 1
    L = compute_Lmod();
else
   L = compute_Lmod(varargin{1});
end

L = double(L);
%save('Ltest','L')

end


function L = compute_Lmod(varargin)
% Compute Lmod matrix for a straight panel of 16 G-L points
%
% This corresponds to determining the weights of I =
% Int(g(t)*log(abs(t-t0))dt,-1,1) (eq. 30 in AK notes)

global high_precision;

% Discretise panel between -1 and 1
[zsrc,worig] = gaussleg(16,[-1 1]); 

% Assign target points
if nargin < 1
    [z1,~] = gaussleg(16,[-3 -1]); %left panel
    [z2,~] = gaussleg(16,[1 3]); %right panel
    ztar = [z1(end-3:end); zsrc; z2(1:4)];
else
    ztar = varargin{1};
end

if high_precision
    L = vpa(zeros(length(ztar),16));
else
    L = zeros(length(ztar),16);
end

for j=1:length(ztar) % Go through all points and compute the weights
    
    % Assign target point
    nz = ztar(j);
    
    if high_precision
        lg1 = vpa(log(abs(1-nz)));
        lg2 = vpa(log(abs(1+nz)));
        % Note the slightly different form here from the standard SQ, this is
        % because we are on the real line, considering only straight panels and
        % on-surface evaluation.
        
        
        % Compute p0 analytically. NB p0 corresponds to p(1)
        p = vpa(zeros(17,1));
        r = vpa(zeros(16,1));
    else
        lg1 = log(abs(1-nz));
        lg2 = log(abs(1+nz));
        p = zeros(17,1);
        r = zeros(16,1);
    end
    
    p(1) = lg1-lg2;
    
    % Compute p and r recursively
    
    for k = 2:17
        % Update p_k
        p(k) = nz*p(k-1) + (1-(-1)^(k-1))/(k-1);
        
        % Update r_{k-1}
        r(k-1) = 1/(k-1)*(lg1-(-1)^(k-1)*lg2-p(k));
        % rk would have an additional scaling, which is not needed as the
        % panel is the canonical panel.
        
    end
    
    % Solve V^T wv = r
    wv = vandernewton(zsrc,r,16); 
    % Also, note here we do not take the imaginary part of r as we're on
    % the real axis.
    
    % Scale the new weights with standard GL weights
    wv = wv./worig;
    
    % Add on additional log from the first sum, log(abs(tl-tj))
    tmp =  ~ismembertol(double(zsrc),double(ztar(j)),1e-10);
    wv(tmp) = wv(tmp) - log(abs(zsrc(tmp)-ztar(j)));
    
    % Add weights to matrix
    L(j,:) = wv.';
end


if nargin < 1
    % Compare with precomputed matrix
    load('LmodMatSame.mat', 'Lmod');
    Lrikard = [Lmod;rot90(Lmod,2)];
    
    fprintf('Difference between L (mine) and Lmod (Rikard): %e \n',double(norm(L-Lrikard,inf)))
end

end

function b = vandernewton(T,b,n)

for k=1:n-1
    for i=n:-1:k+1
        b(i) = b(i) - T(k)*b(i-1);
    end
end

for k=n-1:-1:1
    for i=k+1:n
        b(i) = b(i)/(T(i)-T(i-k));
    end
    for i=k:n-1
        b(i) = b(i) - b(i+1);
    end
end
end


function [nodes,weights] = gaussleg(n,interv)
% Creates G-L nodes and weights for polynomials of degree n, on
% the interval interv (as done in Trefethen)

global high_precision;

n_vec = 1:(n-1);

if high_precision
    interv = vpa(interv);
    beta = vpa(0.5*(1-(2*n_vec).^(-2)).^(-1/2));
else
    beta = 0.5*(1-(2*n_vec).^(-2)).^(-1/2);
end

Tn = diag(beta,-1) + diag(beta,1);
[V,D] = eig(Tn);

if high_precision
    % eigenvalues are sorted and normalized differently in symbolic package apparently
   D = flipud(D(:,end:-1:1));
   V = -V(:,end:-1:1);
end

nodes = D(logical(eye(size(D))));
weights = (2*V(1,:).^2)'; %here we could use our saved W16 weights and
%remap them instead... but we won't save much time doing so

nodes = (interv(1)*(1-nodes)+interv(2)*(1+nodes))/2;
weights =(interv(2)-interv(1))/2*weights;
end