%==========================================================================
% Surrogate-Accelerated Empirical Interpolation for Digital Twin Wave-Based SHM
% Section 5.1: Analytical Case Study – EIM on Double-Gaussian Pulse
% Authors: Abhilash Sreekumar, Linjun Zhong, Dimitrios Chronopoulos
% Date: 2025-04-25
%==========================================================================

close all; clear all; clc;

%% Generate data
nmu         = 29;                              % Number of μ samples
nx          = 400;                             % Grid size per dimension
mu_values   = linspace(0,7,nmu);               % μ ∈ [0,7]
[X, Y]      = meshgrid(linspace(0,20,nx));     % 2D grid
x           = X(:); 
y           = Y(:);
x2plot      = reshape(x,nx,nx); 
y2plot      = reshape(y,nx,nx);

% Preallocate statistics (unused here, placeholder for future metrics)
entropies       = zeros(1,nmu);
variances       = entropies;
std_devs        = entropies;
coeff_variations= entropies;
spatial_autocorrs= entropies;

% Initialize figure for snapshot and residual plots
fig = figure;  
k   = 0;

% Build snapshot matrix
snapshot = zeros(nx^2, nmu);
for imu = 1:nmu
    snapshot(:,imu) = f(x, y, mu_values(imu));
    sol2plot = reshape(snapshot(:,imu), nx, nx);
    
    % Plot only integer μ snapshots
    if mod(mu_values(imu),1)==0
        k = k+1;
        subplot(8,8,k)
        pcolor(x2plot, y2plot, sol2plot);
        shading interp;
        colorbar;
        colormap jet;                     % Use jet colormap
        title(['\mu = ', num2str(mu_values(imu))]);
        axis([0 20 0 20]); 
        axis off; axis square;
    end
    
    % Shannon entropy (not used further here)
    entropies(imu) = -sum(snapshot(:,imu).*log2(snapshot(:,imu)+eps), 'omitnan');
end

% Rank parameters by entropy (high → heterogeneous)
[~, entropy_sorted_indices] = sort(entropies, 'descend');

%% Continuous EIM
iter           = 1;
imustar        = entropy_sorted_indices(1);    % μ* = highest entropy sample
mustar         = mu_values(imustar);
V              = snapshot(:,imustar);           
idx_mp(1)      = find(abs(V)==max(abs(V)),1);  % First interpolation index
idx_mu(1)      = imustar;                      
V(:,1)         = V / V(idx_mp(1));             % Normalize basis w.r.t max

% Build initial collocation matrix
B = V(idx_mp, :);

% Compute initial basis and snapshot errors
for imu = 1:nmu
    F = f(x(idx_mp), y(idx_mp), mu_values(imu));
    c = B \ F;
    err_basis{iter}(:,imu)    = V*c - snapshot(:,imu);
    err_snapshot{iter}(:,imu) = err_basis{iter}(:,imu);
end

% EIM iteration parameters
errTol      = 1e-3;
iterMax     = 50;
maxErr      = inf;
idx_mu_rec  = [];
errorIterPlot = [];

% Greedy EIM loop
while maxErr >= errTol && iter <= iterMax
    disp(['It ',num2str(iter),'/',num2str(iterMax),...
          ', max error = ',num2str(maxErr),...
          ', μ* = μ_',num2str(idx_mu(end)),...
          ' = ',num2str(mu_values(idx_mu(end)))]);
    idx_mu_rec = [idx_mu_rec, idx_mu(end)];
    
    % Plot residual basis fields for integer μ*
    if mod(mu_values(idx_mu(end)),1)==0  
        for imu = 1:nmu
            err2Plot = reshape(err_basis{iter}(:,imu), nx, nx);
            if mod(mu_values(imu),1)==0
                k = k+1;
                subplot(8,8,k)
                pcolor(x2plot, y2plot, err2Plot);
                shading interp; colorbar; colormap jet;
                axis([0 20 0 20]); 
                axis off; axis square;
            end
        end
        errorIterPlot = [errorIterPlot, iter];
    end

    % Find next global max residual (space and parameter)
    [mp_row, mp_col] = find(abs(err_basis{iter})==max(abs(err_basis{iter}(:))),1);
    idx_mp(iter+1) = mp_row;
    idx_mu(iter+1) = mp_col;

    % Enrich basis with normalized selected residual
    newVec = err_basis{iter}(:, idx_mu(end)) / err_basis{iter}(idx_mp(end), idx_mu(end));
    V = [V, newVec];

    % Rebuild collocation matrix
    B = V(idx_mp, :);

    % Update basis errors for each μ
    for imu = 1:nmu
        [found, loc] = ismember(imu, idx_mu);
        if ~found
            F = err_basis{iter}(idx_mp, imu);
            c = B \ F;
            err_basis{iter+1}(:,imu) = V*c - err_basis{iter}(:,imu);
        else
            if loc == 1
                F = V(idx_mp,1);
                c = B \ F;
                err_basis{iter+1}(:,imu) = V*c - V(:,1);
            else
                F = err_basis{loc-1}(idx_mp, imu);
                c = B \ F;
                err_basis{iter+1}(:,imu) = V*c - err_basis{loc-1}(:,imu);
            end
        end
    end

    % Update snapshot errors and max error
    for imu = 1:nmu
        F = snapshot(idx_mp(1:iter+1), imu);
        c = B \ F;
        err_snapshot{iter+1}(:,imu) = V*c - snapshot(:,imu);
    end
    maxErr = max(abs(err_snapshot{iter+1}(:)));
    max_errRecord(iter) = maxErr;

    iter = iter + 1;
end

% Expand figure full-screen
set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

%% Plot EIM convergence & μ-sequence

% Convergence (semi-log y)
figure;
plot(1:length(max_errRecord), max_errRecord, 'b-o','LineWidth',1.5);
set(gca, 'YScale','log','GridLineStyle','--');
grid on;
xlabel('n_b'); ylabel('e^*');
title('EIM – Worst‐case Error vs. # Bases');

% μ_selection trajectory
figure;
plot(1:length(idx_mu_rec), mu_values(idx_mu_rec), 'b-o','LineWidth',1.5);
grid on;
xlabel('n_b'); ylabel('\mu^*');
title('EIM – Selected μ at Each Iteration');

% Spatial scatter of interpolation points colored by μ
figure;
scatter(x(idx_mp), y(idx_mp), 75, mu_values(idx_mu), 'filled');
colormap jet; colorbar; caxis([0,7]);
axis square; grid on;
xlabel('x'); ylabel('y');
title('EIM Interpolation Points (color=μ)');

