function [Pxc, Pyc, Px, Py, X, Y] = evaluate_pressure_gradient(solution, varargin)
%EVALUTATE_PRESSURE evaluates the pressure at either regular grid over the
%reference cell or at specified target points. Uses the stresslet identity
%to identify points inside the domain. Supports periodic and non-periodic 
%domains, and uses special quadrature to evaluate close to the boundary.
%
% inputs:
% -solution: structure containing the following fields
%   -problem: a problem structure containing problem/geometry information
%   -q : density function, vector of size #npts by 2
% variable number of additional arguments, one of
% 1) N: number of points in each direction
% 2) X and Y, vectors (or matrices) of target points

%disp('Evaluating pressure gradient...');

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

if nargin == 2 % N specified, evaluate on regular grid
    N = varargin{1};
    if solution.problem.periodic
        x = linspace(-Lx/2, Lx/2, N);
        y = linspace(-Ly/2, Ly/2, N);
    else
        x = linspace(min(xsrc), max(xsrc), N);
        y = linspace(min(ysrc), max(ysrc), N);
    end
    
    [X,Y] = meshgrid(x,y);
else
    X = varargin{1};
    Y = varargin{2};    
end

if solution.problem.periodic
    pslp = StokesSLP_pressure_grad_ewald_2p(xsrc, ysrc, X(:), Y(:),...
        solution.q(:,1).*weights, solution.q(:,2).*weights, Lx, Ly,...
        'verbose', 0);
    
    % for special quadrature points must be inside reference cell
    Xtar_sq = mod(X+Lx/2,Lx)-Lx/2;
    Ytar_sq = mod(Y+Ly/2,Ly)-Ly/2;
    
    %correct using special quadrature
    [pslp_corrected,~] = mex_SQ_slp_pressure_grad(Xtar_sq(:)+1i*(Ytar_sq(:)+1e-60), ...
        domain.z, domain.zp, domain.quad_weights, domain.panel_breaks,...
        domain.wazp, domain.z32, domain.zp32, domain.quad_weights32,...
        domain.wazp32, solution.q(:,1)+1i*solution.q(:,2),...
        pslp(1,:)'+1i*(pslp(2,:)'+1e-60), domain.mean_panel_length, domain.extra.gridSolidmat, ...
        domain.extra.Nrows, domain.extra.Ncols, domain.extra.panels2wall,...
        domain.reference_cell,true);
    
    if isinf(solution.problem.eta)
        p = pslp;
        p = p(1,:)'+1i*p(2,:)';
        p_corrected = pslp_corrected;
    else
        
        pdlp = StokesDLP_pressure_grad_ewald_2p(xsrc, ysrc, X(:), Y(:), n1, n2,...
            solution.q(:,1).*weights, solution.q(:,2).*weights, Lx, Ly,...
            'verbose', 0);
    
        [pdlp_corrected,~] = mex_SQ_dlp_pressure_grad(Xtar_sq(:)+1i*(Ytar_sq(:)+1e-60),...
            domain.z, domain.zp, domain.quad_weights, domain.panel_breaks,...
            domain.wazp, domain.z32, domain.zp32, domain.quad_weights32,...
            domain.wazp32, solution.q(:,1)+1i*solution.q(:,2), ...
            pdlp(1,:)'+1i*(pdlp(2,:)'+1e-60),...
            domain.mean_panel_length, domain.extra.gridSolidmat, ...
            domain.extra.Nrows, domain.extra.Ncols, domain.extra.panels2wall,...
            domain.reference_cell, true);
        
        p = pdlp + solution.problem.eta*pslp;
        p = p(1,:)'+1i*p(2,:)';
        p_corrected = pdlp_corrected + solution.problem.eta*pslp_corrected;
    end    
else % completion flow
    pdlp = evaluate_double_layer_pressure_gradient_direct(X(:),Y(:),...
        domain.z,domain.zp,solution.q(:,1)+1i*solution.q(:,2),domain.quad_weights);

    %correct using special quadrature
    [pdlp_corrected,~] = mex_SQ_dlp_pressure_grad(X(:)+1i*(Y(:)+1e-60), ...
        domain.z, domain.zp, domain.quad_weights, domain.panel_breaks,...
        domain.wazp, domain.z32, domain.zp32, domain.quad_weights32,...
        domain.wazp32, solution.q(:,1)+1i*solution.q(:,2)+1i*1e-60,...
        real(pdlp)+1i*(imag(pdlp)+1e-60), domain.mean_panel_length, domain.extra.gridSolidmat, ...
        domain.extra.Nrows, domain.extra.Ncols, domain.extra.panels2wall,...
        domain.reference_cell,false);
    
    % add on completion flow gradient pressure Stokeslet
    dPs = completion_contribution_pressure_gradient(domain.centers(2:end),...
        X(:)+1i*Y(:),solution.forces);
    
    p = pdlp + dPs;
    p_corrected = pdlp_corrected + dPs;
end

if solution.problem.periodic && solution.trim
    % find points inside domain by applying stresslet identity
    if solution.problem.periodic
        [test1, test2] = StokesDLP_ewald_2p(xsrc, ysrc, X(:), Y(:), n1, n2,...
            ones(length(n1),1).*weights, zeros(length(n1),1).*weights, Lx, Ly);
    else
        % again, FMM for DLP assumes sources=targets, so we have to artifically
        % add sources with strength zero
        qtmp1 = [ones(length(xsrc),1).*weights; zeros(length(X(:)),1)];
        qtmp2 = [zeros(length(xsrc),1).*weights; zeros(length(X(:)),1)];
        ntmp1 = [n1(:); zeros(length(X(:)),1)];
        ntmp2 = [n2(:); zeros(length(X(:)),1)];

        xtmp = [xsrc; X(:)];
        ytmp = [ysrc; Y(:)];

        [test1, test2] = stokesDLPfmm(qtmp1, qtmp2, xtmp, ytmp, ntmp1, ntmp2);

        % extract data at target points only
        test1 = -test1(length(xsrc)+1:end);
        test2 = -test2(length(xsrc)+1:end);
    end

    % correct using special quadrature
    [test,~] = mex_SQ_dlp(Xtar_sq(:)+1i*(Ytar_sq(:)+1e-60), domain.z,...
                    domain.zp, domain.quad_weights, ...
                    domain.panel_breaks, domain.wazp, domain.z32, domain.zp32,...
                    domain.quad_weights32, domain.wazp32,ones(length(n1),1) + 1e-14*1i,...
                    test1 + 1i*test2,domain.mean_panel_length,domain.extra.gridSolidmat, ...
                    domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
                    domain.reference_cell,solution.problem.periodic);

    % anything that is greater than 0 is outside the fluid domain
    outside = find(solution.problem.stresslet_id_test(real(test)) == 1);
    p_corrected(outside) = nan + 1i*nan;
    p(outside) = nan + 1i*nan;

    X(outside) = nan;
    Y(outside) = nan;
end

dP = reshape(p, size(X));
dPc = reshape(p_corrected, size(X));

Px = real(dP);
Py = imag(dP);
Pxc = real(dPc);
Pyc = imag(dPc);

% add on driving pressure
if solution.problem.periodic
    Px = Px + solution.problem.pressure_gradient_x;
    Py = Py + solution.problem.pressure_gradient_y;
    Pxc = Pxc + solution.problem.pressure_gradient_x;
    Pyc = Pyc + solution.problem.pressure_gradient_y;
end