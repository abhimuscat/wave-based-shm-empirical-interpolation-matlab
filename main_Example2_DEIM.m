%==========================================================================
% Surrogate-Accelerated Empirical Interpolation for Digital Twin Wave-Based SHM
% Section 5.2: DEIM Digital Twin for a Beam with Localized Stiffness Defect
% Authors: Abhilash Sreekumar, Linjun Zhong, Dimitrios Chronopoulos
% Date: 2025-04-25
%==========================================================================
clear all; clc; close all;

%% Training parameter setup
opts.corrFunc    = @corrmatern32;
nmu             = 25;                        % Number of training μ samples
mu_values       = linspace(0,1,nmu);         % μ ∈ [0,1]
nmuTest_local   = 250;                       % # test μ samples (for surrogate)
nb_desired      = 25;                        % Desired number of DEIM bases

%% Simulation parameters
input.L          = 1.0;    % bar length [m]
input.Nx         = 401;    % spatial nodes
input.T          = 0.02;   % total time [s]
input.f0         = 1000;   % central frequency [Hz]
input.rho        = 7800;   % density [kg/m^3]
input.c0         = 300;    % wave speed [m/s]
input.c_d        = 2.5e6;  % damping coeff [N·s/m]
input.damage_width = 0.02; % damage region width [m]
input.D_value     = 0.5;   % 50% stiffness reduction

%% Sensor locations
sensorPositions = [0.2, 0.4, 0.6, 0.8];
Nsensors        = numel(sensorPositions);
x               = linspace(0,input.L,input.Nx)';
sensorIdx       = zeros(Nsensors,1);
for j = 1:Nsensors
    [~, sensorIdx(j)] = min(abs(x - sensorPositions(j)));
end

%% Training simulations
fprintf('Starting training simulations...\n')
tic;
for iMu = 1:nmu
    fprintf(' Training μ-sample %d/%d\n', iMu, nmu);
    input.damage_center = mu_values(iMu);
    tic_local = tic;
    [time, u_record{iMu}] = wavePropagation1d(input);
    time_record(iMu) = toc(tic_local);
end
trainingTime = toc;
fprintf(' Training runtime: %.2f s\n\n', trainingTime);
Nt = numel(time);

%% Extract sensor time series
fprintf('Extracting sensor data...\n')
tic;
for j = 1:Nsensors
    sensorData{j} = zeros(Nt,nmu);
    for iMu = 1:nmu
        sensorData{j}(:,iMu) = u_record{iMu}(sensorIdx(j),:)';
    end
end
extractionTime = toc;
fprintf(' Extraction runtime: %.2f s\n\n', extractionTime);

%% DEIM for each sensor
fprintf('Running DEIM procedure...\n')
tic;
for j = 1:Nsensors
    snapshot = sensorData{j};           % [Nt x nmu]
    
    % POD via SVD
    [Vpod,S,~] = svd(snapshot,'econ');
    Vpod       = Vpod(:,1:nmu);         % take first nmu modes
    
    % Initialize DEIM basis
    V          = Vpod(:,1);
    idx_mp     = find(abs(V)==max(abs(V)),1);
    B          = V(idx_mp,:);
    
    % Initial error computation
    for ib = 1:nmu
        F = Vpod(idx_mp,ib);
        c = B\F;
        err_basis{1}(:,ib)    = V*c - Vpod(:,ib);
        err_snapshot{1}(:,ib) = V*c - snapshot(:,ib);
    end
    
    errTol   = 1e-3;
    iterMax  = nb_desired;
    maxErr   = inf;
    iter     = 1;
    max_errRecord = [];
    
    % Greedy DEIM loop
    while iter < iterMax
        fprintf('  Sensor %d, iter %d, maxErr = %.3e\n', j, iter, maxErr);
        currentErrs = abs(err_basis{iter}(:,iter+1));
        idx_mp_temp = find(currentErrs==max(currentErrs),1);
        idx_mp(end+1) = idx_mp_temp;
        
        V = Vpod(:,1:iter+1);
        B = V(idx_mp,:);
        
        for ib = 1:nmu
            F = Vpod(idx_mp,ib);
            c = B\F;
            err_basis{iter+1}(:,ib) = V*c - Vpod(:,ib);
        end
        for iMu = 1:nmu
            F = snapshot(idx_mp,iMu);
            c = B\F;
            err_snapshot{iter+1}(:,iMu) = V*c - snapshot(:,iMu);
        end
        
        maxErr = max(abs(err_snapshot{iter+1}(:)));
        max_errRecord(iter) = maxErr;
        iter = iter + 1;
    end
    
    % Store DEIM info
    DEIM_result{j}.Vpod = Vpod;
    DEIM_result{j}.idx  = idx_mp;
    max_errRecord_all{j} = max_errRecord;
    
    % Compute DEIM coefficients
    r_DEIM = iter;
    B_final = Vpod(:,1:r_DEIM);
    c_mat   = zeros(r_DEIM,nmu);
    for iMu = 1:nmu
        F = snapshot(idx_mp(1:r_DEIM),iMu);
        c_mat(:,iMu) = B_final(idx_mp(1:r_DEIM),:)\F;
    end
    c_coeff{j} = c_mat;
    
    % Fit Kriging surrogate for each coefficient
    muTest = linspace(0,1,nmuTest_local)';
    cSurr  = zeros(r_DEIM,nmuTest_local);
    for ib = 1:r_DEIM
        model = oodacefit(mu_values.',c_mat(ib,:).',opts);
        cSurr(ib,:) = model.predict(muTest);
    end
    DEIM_surrogate{j}.muTest = muTest;
    DEIM_surrogate{j}.cSurr  = cSurr;
    DEIM_surrogate{j}.Basis  = B_final;
end
deimTime = toc;
fprintf('\nDEIM runtime: %.2f s\n\n', deimTime);

%% Testing with FEM and DEIM surrogate
fprintf('Starting testing simulations...\n')
tic;
muTest = linspace(0,1,nmuTest_local);
for iMu = 1:nmuTest_local
    fprintf(' Test μ-sample %d/%d\n', iMu, nmuTest_local);
    input.damage_center = muTest(iMu);
    [~,u_test] = wavePropagation1d(input);
    for j = 1:Nsensors
        FEMsig = u_test(sensorIdx(j),:)';
        snapshot_FEM{j}(:,iMu) = FEMsig;
        c_hat = DEIM_surrogate{j}.cSurr(:,iMu);
        SURRsig = DEIM_surrogate{j}.Basis * c_hat;
        snapshot_SURR{j}(:,iMu) = SURRsig;
    end
end
testingTime = toc;
fprintf('\nTesting runtime: %.2f s\n\n', testingTime);

%% Global error metrics
globalErrSum = 0; count = 0;
for j = 1:Nsensors
    relErrs = vecnorm(snapshot_FEM{j}-snapshot_SURR{j}) ./ vecnorm(snapshot_FEM{j});
    globalErrSum = globalErrSum + sum(relErrs);
    count = count + numel(relErrs);
    errors(j,:) = relErrs;
end
avgError = globalErrSum/count;
maxError = max(errors(:));
fprintf('Global avg rel error: %.3e\n', avgError);
fprintf('Global max rel error: %.3e\n', maxError);

%% Plotting – single sensor and μ
selSensor = 3;
selMuIdx   = round(nmuTest_local/3);
figure;
plot(time, snapshot_FEM{selSensor}(:,selMuIdx),'b-','LineWidth',2); hold on;
plot(time, snapshot_SURR{selSensor}(:,selMuIdx),'r--','LineWidth',2);
xlabel('Time [s]'); ylabel('Displacement');
title(sprintf('Sensor x=%.2f, μ=%.2f',sensorPositions(selSensor),muTest(selMuIdx)));
legend('FEM','DEIM','Location','Best');

%% Plot DEIM convergence per sensor
figure;
for j = 1:Nsensors
    subplot(2,ceil(Nsensors/2),j);
    semilogy(max_errRecord_all{j}(1:end-1),'LineWidth',2);
    xlabel('Iteration'); ylabel('Max Error');
    title(sprintf('Sensor x=%.2f',sensorPositions(j)));
    grid on;
end
sgtitle('DEIM Convergence');

%% Plot surrogate vs training coefficients for one sensor
nb = size(c_coeff{selSensor},1);
nrows = floor(sqrt(nb)); ncols = ceil(nb/nrows);
figure;
for ib = 1:nb
    subplot(nrows,ncols,ib);
    plot(DEIM_surrogate{selSensor}.muTest, DEIM_surrogate{selSensor}.cSurr(ib,:), 'r-','LineWidth',1.5); hold on;
    plot(mu_values, c_coeff{selSensor}(ib,:),'bo','MarkerSize',5);
    xlabel('\mu'); ylabel(sprintf('c_%d',ib));
    title(sprintf('Coeff %d',ib));
    legend('Surrogate','Training','Location','Best');
    grid on;
end
sgtitle(sprintf('Coefficients for Sensor x=%.2f',sensorPositions(selSensor)));


