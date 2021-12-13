function [omegac, omega, X, Y] = evaluate_vorticity(solution, varargin)

domain = solution.problem.domain;
if solution.problem.periodic
    Lx = domain.Lx;
    Ly = domain.Ly;
end

xsrc = real(solution.problem.domain.z);
ysrc = imag(solution.problem.domain.z);
n1 = real(-1i*solution.problem.domain.zp)./abs(solution.problem.domain.zp);
n2 = imag(-1i*solution.problem.domain.zp)./abs(solution.problem.domain.zp);
weights = solution.problem.domain.wazp;

% default parameters
fmm = 1;
verbose = 0;

if nargin == 2 || (nargin > 2 && ischar(varargin{2})) % N specified, evaluate on regular grid
    N = varargin{1};
    if solution.problem.periodic
        x = linspace(-Lx/2, Lx/2, N);
        y = linspace(-Ly/2, Ly/2, N);
    else
        x = linspace(min(xsrc), max(xsrc), N);
        y = linspace(min(ysrc), max(ysrc), N);
    end
    
    [X,Y] = meshgrid(x,y);
elseif nargin == 3 && ~ischar(varargin{2}) % target points are specified
    X = varargin{1};
    Y = varargin{2};
elseif nargin > 3 
    % given target points
    X = varargin{1};
    Y = varargin{2};

    % Go through all other input arguments and assign parameters
    jv = 3;
    while jv <= length(varargin)-1
       switch varargin{jv}
           case 'fmm'
               fmm = varargin{jv+1}; 
           case 'verbose'
               verbose = varargin{jv+1};
       end
       jv = jv + 2;
    end
end

if solution.problem.periodic
    omega_slp = StokesSLP_vorticity_ewald_2p(xsrc, ysrc, X(:), Y(:),...
        solution.q(:,1).*weights, solution.q(:,2).*weights, Lx, Ly,...
        'verbose', verbose);
    
    % for special quadrature points must be inside reference cell
    Xtar_sq = mod(X+Lx/2,Lx)-Lx/2;
    Ytar_sq = mod(Y+Ly/2,Ly)-Ly/2;
    
    [omega_slp_corrected,~] = mex_SQ_slp_vorticity(Xtar_sq(:)+1i*(Ytar_sq(:)+1e-60), ...
        domain.z, domain.zp, domain.quad_weights, domain.panel_breaks,...
        domain.wazp, domain.z32, domain.zp32, domain.quad_weights32,...
        domain.wazp32, solution.q(:,1)+1i*solution.q(:,2),...
        omega_slp, domain.mean_panel_length, domain.extra.gridSolidmat, ...
        domain.extra.Nrows, domain.extra.Ncols, domain.extra.panels2wall,...
        domain.reference_cell,true);
    
    if isinf(solution.problem.eta)
        omega = omega_slp;
        omega_corrected = omega_slp_corrected;
    else
        omega_dlp = StokesDLP_vorticity_ewald_2p(xsrc, ysrc, X(:), Y(:), n1, n2,...
            solution.q(:,1).*weights, solution.q(:,2).*weights, Lx, Ly,...
            'verbose', verbose);
        
        [omega_dlp_corrected,~] = mex_SQ_dlp_vorticity(Xtar_sq(:)+1i*(Ytar_sq(:)+1e-60), ...
            domain.z, domain.zp, domain.quad_weights, domain.panel_breaks, domain.wazp, domain.z32,...
            domain.zp32, domain.quad_weights32, domain.wazp32,...
            solution.q(:,1)+1i*solution.q(:,2),omega_dlp,domain.mean_panel_length,...
            domain.extra.gridSolidmat,domain.extra.Nrows,domain.extra.Ncols,...
            domain.extra.panels2wall,domain.reference_cell, true);
        
        omega = omega_dlp + solution.problem.eta*omega_slp;
        omega_corrected = omega_dlp_corrected + solution.problem.eta*omega_slp_corrected;
    end
else
    % if the problem is not periodic, use FMM and the completed double-layer
    % potential, i.e. Power-Miranda (dlp + stokeslet + rotlet)
    omega_dlp = evaluate_double_layer_vorticity_direct(X(:),Y(:),...
        domain.z,domain.zp,solution.q(:,1)+1i*solution.q(:,2),domain.quad_weights);
    
    % correct using special quadrature
    [omega_dlp_corrected,~] = mex_SQ_dlp_vorticity(X(:)+1i*(Y(:)+1e-60), domain.z, domain.zp,...
        domain.quad_weights, domain.panel_breaks, domain.wazp, domain.z32,...
        domain.zp32, domain.quad_weights32, domain.wazp32,...
        solution.q(:,1)+1i*solution.q(:,2),real(omega_dlp)+1i*(imag(omega_dlp)+1e-60),domain.mean_panel_length,...
        domain.extra.gridSolidmat,domain.extra.Nrows,domain.extra.Ncols,...
        domain.extra.panels2wall,domain.reference_cell, false);
    
    % add on completion flow from rotlets and Stokeslets
    omegaR = zeros(length(X(:)),1);
    omegaS = completion_contribution_vorticity(domain.centers(2:end),X(:)+1i*Y(:),...
        solution.forces);
    
    omega = omega_dlp + omegaS + omegaR;
    omega_corrected = omega_dlp_corrected + omegaS + omegaR;
end

omega = real(omega);
omega_corrected = real(omega_corrected);

if solution.problem.periodic && solution.trim
    % find points inside domain by applying stresslet identity
end

omega = reshape(omega, size(X));
omegac = reshape(omega_corrected, size(X));