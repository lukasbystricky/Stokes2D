function [U1c, U2c, X, Y, U1, U2] = evaluate_velocity(solution, varargin)
%EVALUTATE_VELOCITY evaluates the velocity at either regular grid over the 
%reference cell or at specified target points. Uses the stresslet identity 
%to identify points inside the domain.
%
% inputs:
% -solution: structure containing the following fields
%   -problem: a problem structure containing problem/geometry information
%   -q : density function, vector of size #npts by 2
% variable number of additional arguments, one of
% 1) N: number of points in each direction
% 2) X and Y, vectors (or matrices) of target points

disp('Evaluating velocity...');

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
else % target points are specified
    X = varargin{1};
    Y = varargin{2};
end

% if the problem is periodic use spectral Ewald and combined layer
% formulation
if solution.problem.periodic
   
    [uslp1, uslp2] = StokesSLP_ewald_2p(xsrc, ysrc, X(:), Y(:),...
        solution.q(:,1).*weights, solution.q(:,2).*weights, Lx, Ly,...
        'verbose', 1);
    
    uslp = uslp1 + 1i*uslp2;
    if ~isinf(solution.problem.eta)
        
        [udlp1, udlp2] = StokesDLP_ewald_2p(xsrc, ysrc, X(:), Y(:), n1, n2,...
            solution.q(:,1).*weights, solution.q(:,2).*weights, Lx, Ly,...
            'verbose', 1);
        udlp = udlp1 + 1i*udlp2;
        
    end
        
    % for special quadrature points must be inside reference cell
    Xtar_sq = mod(X+Lx/2,Lx)-Lx/2;
    Ytar_sq = mod(Y+Ly/2,Ly)-Ly/2;
    Xsrc_sq = real(domain.z);
    Ysrc_sq = imag(domain.z);
    
    % correct using special quadrature
    [uslp_corrected,~] = mex_SQ_slp(Xtar_sq(:)+1i*Ytar_sq(:), Xsrc_sq+1i*Ysrc_sq,...
        domain.zp, domain.quad_weights, domain.panel_breaks, domain.wazp, domain.z32,...
        domain.zp32, domain.quad_weights32, domain.wazp32, ...
        solution.q(:,1)+1i*solution.q(:,2),...
        uslp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
        domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
        domain.reference_cell, true);
    
    if isinf(solution.problem.eta)
        u_corrected = uslp_corrected  + ...
            solution.u_avg(1) + 1i*solution.u_avg(2);
        u = uslp + solution.u_avg(1) + 1i*solution.u_avg(2);
    else
        [udlp_corrected,~] = mex_SQ_dlp(X(:)+1i*Y(:), domain.z, domain.zp,...
            domain.quad_weights, domain.panel_breaks, domain.wazp, domain.z32,...
            domain.zp32, domain.quad_weights32, domain.wazp32,...
            solution.q(:,1)+1i*solution.q(:,2),...
            udlp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
            domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
            domain.reference_cell, true);
        
        u_corrected = udlp_corrected + solution.problem.eta*uslp_corrected  + ...
            solution.u_avg(1) + 1i*solution.u_avg(2);
        u = udlp + solution.problem.eta*uslp + solution.u_avg(1) +...
            1i*solution.u_avg(2);
    
    end
        
% if the problem is not periodic, use FMM and the completed double-layer
% potential, i.e. Power-Miranda
else
    
    % FMM can only evaluate for source=targets, so we include all the
    % target points as source points with strength 0
    
    qtmp1 = [solution.q(:,1).*weights; zeros(length(X(:)),1)];
    qtmp2 = [solution.q(:,2).*weights; zeros(length(X(:)),1)];
    ntmp1 = [n1(:); zeros(length(X(:)),1)];
    ntmp2 = [n2(:); zeros(length(X(:)),1)];
    
    xtmp = [xsrc; X(:)];
    ytmp = [ysrc; Y(:)];
    
    [udlp1, udlp2] = stokesDLPfmm(qtmp1, qtmp2, xtmp, ytmp, ntmp1, ntmp2);
    
    % note negative sign in front of double-layer
    udlp = -udlp1(length(xsrc)+1:end) - 1i*udlp2(length(xsrc)+1:end);
    
    disp('Beginning special quadrature...');
    
    % correct using special quadrature
    [udlp_corrected,~] = mex_SQ_dlp(X(:)+1i*Y(:), domain.z, domain.zp,...
        domain.quad_weights, domain.panel_breaks, domain.wazp, domain.z32,...
        domain.zp32, domain.quad_weights32, domain.wazp32,...
        solution.q(:,1)+1i*solution.q(:,2),...
        udlp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
        domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
        domain.reference_cell, false);
    
    % add on completion flow from rotlets and Stokeslets
    [uS, uR] = completion_contribution(domain.centers(2:end), X(:)+1i*Y(:),...
                solution.forces, solution.torques);
    
    u_corrected = udlp_corrected + uS + uR;
    u = udlp + uS + uR;

end

u1 = real(u);
u2 = imag(u);
u1_corrected = real(u_corrected);
u2_corrected = imag(u_corrected);

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
[test,~] = mex_SQ_dlp(X(:)+1i*Y(:), domain.z, domain.zp, domain.quad_weights, ...
                domain.panel_breaks, domain.wazp, domain.z32, domain.zp32,...
                domain.quad_weights32, domain.wazp32,ones(length(n1),1) + 1e-14*1i,...
                test1 + 1i*test2,domain.mean_panel_length,domain.extra.gridSolidmat, ...
                domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
                domain.reference_cell,solution.problem.periodic);

% anything that is greater than 0 is outside the fluid domain
outside = find(solution.problem.stresslet_id_test(real(test)) == 1);
u1_corrected(outside) = nan;
u2_corrected(outside) = nan;

u1(outside) = nan;
u2(outside) = nan;

% X(outside) = nan;
% Y(outside) = nan;

U1c = reshape(u1_corrected, size(X));
U2c = reshape(u2_corrected, size(X));

U1 = reshape(u1, size(X));
U2 = reshape(u2, size(X));