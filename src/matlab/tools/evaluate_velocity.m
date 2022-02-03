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

%disp('Evaluating velocity...');

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

% if the problem is periodic use spectral Ewald and combined layer
% formulation
if solution.problem.periodic
   
    [uslp1, uslp2] = StokesSLP_ewald_2p(xsrc, ysrc, X(:), Y(:),...
        solution.q(:,1).*weights, solution.q(:,2).*weights, Lx, Ly,...
        'verbose', verbose);
    
    uslp = uslp1 + 1i*uslp2;
    if ~isinf(solution.problem.eta)
        
        [udlp1, udlp2] = StokesDLP_ewald_2p(xsrc, ysrc, X(:), Y(:), n1, n2,...
            solution.q(:,1).*weights, solution.q(:,2).*weights, Lx, Ly,...
            'verbose', verbose);
        udlp = udlp1 + 1i*udlp2;
        
    end
        
    % for special quadrature points must be inside reference cell
    Xtar_sq = mod(X+Lx/2,Lx)-Lx/2;
    Ytar_sq = mod(Y+Ly/2,Ly)-Ly/2;
    Xsrc_sq = real(domain.z);
    Ysrc_sq = imag(domain.z);
    
    % correct using special quadrature
    [uslp_corrected,~] = mex_SQ_slp(Xtar_sq(:)+1i*(Ytar_sq(:)+1e-60), Xsrc_sq+1i*(Ysrc_sq+1e-60),...
        domain.zp, domain.quad_weights, domain.panel_breaks, domain.wazp, domain.z32,...
        domain.zp32, domain.quad_weights32, domain.wazp32, ...
        solution.q(:,1)+1i*(solution.q(:,2)+1e-60),...
        uslp,domain.mean_panel_length,domain.extra.gridSolidmat, ...
        domain.extra.Nrows,domain.extra.Ncols,domain.extra.panels2wall,...
        domain.reference_cell, true);
    
    if isinf(solution.problem.eta)
        u_corrected = uslp_corrected  + ...
            solution.u_avg(1) + 1i*solution.u_avg(2);
        u = uslp + solution.u_avg(1) + 1i*solution.u_avg(2);
    else
        [udlp_corrected,~] = mex_SQ_dlp(Xtar_sq(:)+1i*(Ytar_sq(:)+1e-60), Xsrc_sq+1i*(Ysrc_sq+1e-60),...
            domain.zp, domain.quad_weights, domain.panel_breaks, domain.wazp, domain.z32,...
            domain.zp32, domain.quad_weights32, domain.wazp32,...
            solution.q(:,1)+1i*(solution.q(:,2)+1e-60),...
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
    if fmm
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
    
    else
        udlp1 = zeros(numel(X(:)),1);
        udlp2 = zeros(numel(X(:)),1);
        for k = 1:numel(X(:))
            rx = X(k) - xsrc;
            ry = Y(k) - ysrc;
            rho4 = (rx.^2 + ry.^2).^2;
            
            rdotq = rx.*solution.q(:,1).*weights + ry.*solution.q(:,2).*weights;
            rdotn = rx.*n1 + ry.*n2;
            
            udlp1(k) = 4*sum(rdotn.*rdotq./rho4.*rx);
            udlp2(k) = 4*sum(rdotn.*rdotq./rho4.*ry);
        end
        
        udlp1 = udlp1/4/pi;
        udlp2 = udlp2/4/pi;
        
        udlp = -udlp1 - 1i*udlp2;
    end
    disp('Beginning special quadrature...');
    
    % correct using special quadrature
    [udlp_corrected,~] = mex_SQ_dlp(X(:)+1i*(Y(:)+1e-60), domain.z, domain.zp,...
        domain.quad_weights, domain.panel_breaks, domain.wazp, domain.z32,...
        domain.zp32, domain.quad_weights32, domain.wazp32,...
        solution.q(:,1)+1i*(solution.q(:,2)+1e-60),...
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

if solution.problem.periodic && solution.trim
    % find points inside domain by applying stresslet identity
    if solution.problem.periodic
        [test1, test2] = StokesDLP_ewald_2p(xsrc, ysrc, X(:), Y(:), n1, n2,...
            ones(length(n1),1).*weights, zeros(length(n1),1).*weights, Lx, Ly);
    else
        if fmm
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
        else
            test1 = zeros(numel(q1),1);
            test2 = zeros(numel(q1),1);
            for k = 1:numel(q1)


                rx = xsrc(k) - xsrc;
                ry = ysrc(k) - ysrc;
                rho4 = (rx.^2 + ry.^2).^2;

                rdotq = rx.*ones(numel(q1)).*weights;
                rdotn = rx.*n1 + ry.*n2;

                test1(k) = 4*sum(rdotn.*rdotq./rho4.*rx);
                test2(k) = 4*sum(rdotn.*rdotq./rho4.*ry);
            end

            test1 = test1/4/pi;
            test2 = test2/4/pi;         
        end
    end

    %correct using special quadrature
    [test,~] = mex_SQ_dlp(Xtar_sq(:)+1i*(Ytar_sq(:)+1e-60), Xsrc_sq+1i*(Ysrc_sq+1e-60),...
                    domain.zp, domain.quad_weights, ...
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

    X(outside) = nan;
    Y(outside) = nan;
end

U1c = reshape(u1_corrected, size(X));
U2c = reshape(u2_corrected, size(X));

U1 = reshape(u1, size(X));
U2 = reshape(u2, size(X));