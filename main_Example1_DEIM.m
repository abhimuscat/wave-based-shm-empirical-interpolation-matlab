%==========================================================================
% Surrogate-Accelerated Empirical Interpolation for Digital Twin Wave-Based SHM
% Section 5.1: Analytical Case Study – DEIM on Double-Gaussian Pulse
% Authors: Abhilash Sreekumar, Linjun Zhong, Dimitrios Chronopoulos
% Date: 2025-04-25
%==========================================================================

% close all;          % Uncomment to close all figures at start
clear all; clc;       % Clear workspace and command window

%% Generate data
nmu = 29;                                  % Number of parameter samples (μ)
nx  = 400;                                 % Number of spatial grid points per dimension
mu_values = linspace(0,7,nmu);             % μ ∈ [0,7]
[X, Y] = meshgrid(linspace(0,20,nx));      % 2D spatial grid on [0,20]×[0,20]
x = X(:);  y = Y(:);                       % Flatten into vectors
x2plot = reshape(x,nx,nx); 
y2plot = reshape(y,nx,nx);

% Preallocate snapshot matrix
snapshot = zeros(nx^2, nmu);

% Initialize plot counter and figure
fig = figure;  
k = 0;

% Generate and plot snapshots for integer μ values
for imu = 1:nmu
    % Evaluate the double-Gaussian field at this μ
    snapshot(:,imu) = f(x, y, mu_values(imu));
    
    % Reshape for plotting
    sol2plot = reshape(snapshot(:,imu), nx, nx);

    if mod(mu_values(imu),1)==0      % Only if μ is integer
        k = k + 1;
        subplot(8,8,k)
        pcolor(x2plot, y2plot, sol2plot);
        shading interp;
        colorbar;
        xlabel('x'); ylabel('y');
        title(['\mu = ', num2str(mu_values(imu))]);
        axis([0 20 0 20]);
        axis off; axis square;
    end
end

%% DEIM Reduction
% 1) POD via SVD ("econ" for economy size)
[Vpod, SIG, ~] = svd(snapshot, "econ");
SIG = diag(SIG);
nb  = nmu;               % Number of POD modes to keep
Vpod = Vpod(:,1:nb);      % Truncate POD basis

% 2) Initialize DEIM indices and error storage
iter      = 1;
V         = Vpod(:,1);
idx_mp(1) = find(abs(V)==max(abs(V)),1);  % First DEIM index
B         = V(idx_mp, :);                 % Collocation matrix

% Compute initial basis and snapshot residuals
for ib = 1:nb
    F = Vpod(idx_mp, ib);
    c = B \ F;
    err_basis{iter}(:,ib) = V*c - Vpod(:,ib);
end

for imu = 1:nmu
    F = snapshot(idx_mp, imu);
    c = B \ F;
    err_snapshot{iter}(:,imu) = V*c - snapshot(:,imu);
end

% DEIM iteration parameters
iterList     = [1,6,10,15,16,17,25];
errTol       = 1e-4;
iterMax      = 50;
maxErr       = inf;

% 3) DEIM greedy loop
while maxErr >= errTol && iter <= iterMax
    disp(['It ',num2str(iter),'/',num2str(iterMax),...
          ' — max. error = ',num2str(maxErr)]);
    
    % Plot residual snapshots at selected iterations
    if ismember(iter, iterList)
        for imu = 1:nmu
            if mod(mu_values(imu),1)==0
                k = k + 1;
                subplot(8,8,k)
                err2Plot = reshape(err_snapshot{iter}(:,imu), nx, nx);
                pcolor(x2plot, y2plot, err2Plot);
                shading interp;
                colorbar; colormap jet;
                axis([0 20 0 20]);
                axis off; axis square;
            end
        end
    end

    % Select next DEIM index from the next basis residual
    r = err_basis{iter}(:, iter+1);
    idx_mp(iter+1) = find(abs(r)==max(abs(r)),1);
    
    % Update collocation basis
    V = Vpod(:,1:iter+1);
    B = V(idx_mp, :);
    
    % Recompute basis residuals
    for ib = 1:nb
        F = Vpod(idx_mp, ib);
        c = B \ F;
        err_basis{iter+1}(:,ib) = V*c - Vpod(:,ib);
    end
    
    % Recompute snapshot residuals
    for imu = 1:nmu
        F = snapshot(idx_mp, imu);
        c = B \ F;
        err_snapshot{iter+1}(:,imu) = V*c - snapshot(:,imu);
    end
    
    % Update maximum error and record
    maxErr = max(abs(err_snapshot{iter+1}(:)));
    max_errRecord(iter) = maxErr;
    iter = iter + 1;
end

% Expand the initial figure to full screen
set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

%% Plot DEIM convergence and selected points

% Convergence curve (log-log)
figure;
plot(1:iter-2, max_errRecord(1:end-1), 'r-o', 'LineWidth', 1.5);
set(gca, 'XScale','lin', 'YScale','log', 'GridLineStyle','--');
grid on;
title('DEIM – Maximum Error vs. Number of Bases');
xlabel('# bases');
ylabel('max error');

% Scatter of DEIM-selected interpolation points
figure;
plot(x(idx_mp), y(idx_mp), 'ko', 'LineWidth', 1.5);
axis([0 20 0 20]); axis square; grid on;
xlabel('x');
ylabel('y');
title('DEIM Selected Interpolation Points');
