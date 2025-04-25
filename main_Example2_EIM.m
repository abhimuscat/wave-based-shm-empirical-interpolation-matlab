%==========================================================================
% Surrogate-Accelerated Empirical Interpolation for Digital Twin Wave-Based SHM
% Section 5.2: EIM Digital Twin for a Beam with Localized Stiffness Defect
% Authors: Abhilash Sreekumar, Linjun Zhong, Dimitrios Chronopoulos
% Date: 2025-04-25
%==========================================================================
clear all; clc; close all;

%% 1. Training Setup
opts.corrFunc     = @corrmatern32;    % Kriging correlation function
nmu               = 25;               % Number of training μ samples
mu0               = 0; mu1 = 1;
mu_values         = linspace(mu0, mu1, nmu);
nmuTest_local     = 250;              % Number of test μ samples
nb_desired        = 25;               % Desired reduced basis size (≤ nmu)

%% 2. Simulation Parameters
input.L           = 1.0;              % Bar length [m]
input.Nx          = 401;              % Spatial nodes
x                 = linspace(0,input.L,input.Nx)';

input.T           = 0.02;             % Simulation time [s]
input.f0          = 1e3;              % Ricker central frequency [Hz]
input.rho         = 7800;             % Density [kg/m^3]
input.c0          = 300;              % Undamaged wave speed [m/s]
input.c_d         = 2.5e6;            % Damping coeff [N·s/m]
input.damage_width= 0.02;             % Damage region width [m]
input.D_value     = 0.5;              % 50% stiffness reduction

%% 3. Sensor Locations
sensorPositions   = [0.2,0.4,0.6,0.8];
Nsensors          = numel(sensorPositions);
sensorIdx         = arrayfun(@(s) find(abs(x-s)==min(abs(x-s)),1), sensorPositions)';

%% 4. Generate Training Snapshots
fprintf('► Running %d training simulations...\n', nmu);
tic;
u_record = cell(1,nmu);
for i = 1:nmu
    fprintf('  – μ = %.3f (%d/%d)\n', mu_values(i), i, nmu);
    input.damage_center = mu_values(i);
    [~, u_record{i}]    = wavePropagation1d(input);
end
trainingTime = toc;
fprintf('  → Training time: %.2f s\n\n', trainingTime);

%% 5. Extract Sensor Time-Series
fprintf('► Extracting sensor data...\n');
tic;
sensorData = cell(1,Nsensors);
Nt = size(u_record{1},2);
for j = 1:Nsensors
    M = zeros(Nt,nmu);
    for i = 1:nmu
        M(:,i) = u_record{i}(sensorIdx(j),:)';
    end
    sensorData{j} = M;
end
extractionTime = toc;
fprintf('  → Extraction time: %.2f s\n\n', extractionTime);

%% 6. EIM Basis Construction & Surrogate Training
fprintf('► Building EIM surrogates (each sensor)...\n');
tic;
Vcell        = cell(1,Nsensors);
Bcell        = cell(1,Nsensors);
idxcell      = cell(1,Nsensors);
errcell      = cell(1,Nsensors);
maxErrRec    = cell(1,Nsensors);
surrogate    = cell(1,Nsensors);
coeff_train  = cell(1,Nsensors);

for j = 1:Nsensors
    S       = sensorData{j};          % [Nt x nmu]
    imustar = ceil(nmu/2);            % reference snapshot index

    % Initialize basis and pivots
    V       = S(:,imustar);
    [~, p]  = max(abs(V));             
    idxs    = p;
    V       = V / V(p);
    B       = V(p,:);

    % Compute initial error
    err    = zeros(size(S));
    for i = 1:nmu
        c = B \ S(p,i);
        err(:,i) = V*c - S(:,i);
    end

    % Greedy EIM loop
    maxErr = max(abs(err(:)));
    rec    = maxErr;
    iter   = 1;
    while iter < nb_desired
        fprintf('  Sensor %d: iter %d, maxErr = %.3e\n', j, iter, maxErr);

        [~, idxLin]       = max(abs(err(:)));
        [pNew, colNew]    = ind2sub(size(err), idxLin);
        idxs(iter+1)      = pNew;

        vNew              = err(:,colNew);
        vNew              = vNew / vNew(pNew);
        V                 = [V, vNew];
        B                 = V(idxs,:);

        % update error
        errNew = zeros(size(err));
        for i = 1:nmu
            c = B \ S(idxs,i);
            errNew(:,i) = V*c - S(:,i);
        end
        err    = errNew;
        maxErr = max(abs(err(:)));
        rec(end+1)= maxErr;
        iter  = iter + 1;
    end

    % Store convergence
    maxErrRec{j} = rec;

    % Final coefficients
    r           = size(V,2);
    Ctrain      = zeros(r,nmu);
    for i = 1:nmu
        Ctrain(:,i) = B \ S(idxs,i);
    end
    coeff_train{j} = Ctrain;

    % Train Kriging surrogates for each coefficient
    muTest     = linspace(mu0,mu1,nmuTest_local)';
    Csur       = zeros(r,nmuTest_local);
    for k = 1:r
        model     = oodacefit(mu_values.', Ctrain(k,:).', opts);
        Csur(k,:) = model.predict(muTest);
    end

    % Save results
    Vcell{j}     = V;
    Bcell{j}     = B;
    idxcell{j}   = idxs;
    surrogate{j}.muTest = muTest;
    surrogate{j}.Csur   = Csur;
end
eimTime = toc;
fprintf('  → EIM & surrogate time: %.2f s\n\n', eimTime);

%% 7. Testing & Reconstruction
fprintf('► Running %d testing simulations...\n', nmuTest_local);
tic;
snapshot_FE = cell(1,Nsensors);
snapshot_SU = cell(1,Nsensors);

for j = 1:Nsensors
    snapshot_FE{j} = zeros(Nt,nmuTest_local);
    snapshot_SU{j} = zeros(Nt,nmuTest_local);
end

muTest = linspace(mu0,mu1,nmuTest_local);
for i = 1:nmuTest_local
    input.damage_center = muTest(i);
    [~, utest]          = wavePropagation1d(input);
    for j = 1:Nsensors
        Sfull                = utest(sensorIdx(j),:)';
        snapshot_FE{j}(:,i)  = Sfull;
        cHat                  = surrogate{j}.Csur(:,i);
        snapshot_SU{j}(:,i)   = Vcell{j} * cHat;
    end
end
testTime = toc;
fprintf('  → Testing time: %.2f s\n\n', testTime);

%% 8. Global Error Metrics
allErr = zeros(Nsensors,nmuTest_local);
for j = 1:Nsensors
    for i = 1:nmuTest_local
        e    = snapshot_FE{j}(:,i)-snapshot_SU{j}(:,i);
        allErr(j,i) = norm(e)/norm(snapshot_FE{j}(:,i));
    end
end
avgErr = mean(allErr(:));
maxErr = max(allErr(:));
fprintf('Global avg rel error: %.3e\n',avgErr);
fprintf('Global max rel error: %.3e\n\n',maxErr);

%% 9. Plots
% (a) One sensor & μ-time trace
selJ = 3; selI = ceil(nmuTest_local/3);
figure; hold on;
plot((1:Nt)*(input.T/Nt), snapshot_FE{selJ}(:,selI),'b-','LineWidth',1.5);
plot((1:Nt)*(input.T/Nt), snapshot_SU{selJ}(:,selI),'r--','LineWidth',1.5);
xlabel('Time [s]'), ylabel('Disp.'); 
title(sprintf('Sensor %.2f m, μ=%.3f', sensorPositions(selJ), muTest(selI)));
legend('Full FE','EIM Surrogate');

% (b) EIM convergence per sensor
figure;
for j=1:Nsensors
    subplot(2,2,j);
    semilogy(maxErrRec{j}(1:end-1),'LineWidth',1.2);
    xlabel('Iteration'), ylabel('Max Err');
    title(sprintf('Sensor @ %.2f m', sensorPositions(j)));
end
sgtitle('EIM Max-Error Convergence');

% (c) Coefficient fits for selJ
r = size(coeff_train{selJ},1);
nrows = floor(sqrt(r)); ncols = ceil(r/nrows);
figure;
for k=1:r
    subplot(nrows,ncols,k); hold on;
    plot(surrogate{selJ}.muTest, surrogate{selJ}.Csur(k,:),'r-');
    plot(mu_values, coeff_train{selJ}(k,:), 'bo');
    xlabel('\mu'), ylabel(sprintf('c_{%d}',k));
    title(sprintf('Coef %d',k));
    legend('Surrogate','Train','Location','Best');
end
sgtitle(sprintf('EIM Coefficients @ Sensor %.2f m', sensorPositions(selJ)));

