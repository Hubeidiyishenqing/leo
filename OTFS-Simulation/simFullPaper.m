%--------------------------------------------------------------------------
%
%  Comprehensive Simulation for Paper:
%    "LEO-NTN OTFS Integrated Communication and Navigation System Design"
%
%  Three Innovation Points:
%    1. OTFS Integrated Communication & Navigation (ISAC) via DD Pilot
%    2. DD-Domain Pilot + Doppler Pre-compensation Joint Design
%    3. OFDM / OTFS / OTFS-Pilot Three-Scheme LEO Performance Comparison
%
%  Generates 4 Figures:
%    Fig 1: Three-scheme BER vs Eb/No (baseline)
%    Fig 4: Pilot Overhead vs Pre-compensation Ratio (analytical)
%    Fig 5: Velocity/Range RMSE vs SNR (with CRB)
%    Fig 7: DD Grid Pilot Pattern Visualization
%
%--------------------------------------------------------------------------

clear; close all;
rng(42);
tTotal = tic;

% Initialize parallel pool (if available)
if isempty(gcp('nocreate'))
    try
        parpool;
        fprintf('Parallel pool started with %d workers.\n', gcp().NumWorkers);
    catch
        fprintf('Parallel Computing Toolbox not available. Running sequentially.\n');
    end
end

%% ========================================================================
%  System Parameters (consistent with main.m)
%  ========================================================================
ModOrder = 16;                      % 16-QAM
k = log2(ModOrder);
cpSize = 0.07;
scs = 120e3;                        % 120 kHz SCS
Bw = 100e6;
ofdmSym = 128;                      % Doppler bins M
velocity = 26000;                    % km/hr (LEO)
fc = 2e9;                           % S-band
codeRate = 1/2;                     % LDPC code rate
maxLDPCIter = 25;                   % LDPC decoder max iterations
c_light = physconst('LightSpeed');

numSC = pow2(ceil(log2(Bw/scs)));    % Delay bins N = 1024
cpLen = floor(cpSize * numSC);
numDC = numSC - 12;                  % Data carriers = 1012

% DD grid resolutions
v_ms = velocity * 1e3 / 3600;
fd_hz = v_ms * fc / c_light;
Ts_sym = (1 + cpSize) / scs;
delta_tau_s = 1 / (numSC * scs);
delta_nu_hz = 1 / (ofdmSym * Ts_sym);

% Default NTN config
ntnConfig.profile     = 'NTN-TDL-C';
ntnConfig.elevAngle   = 50;
ntnConfig.scenario    = 'Suburban';
ntnConfig.fc          = fc;
ntnConfig.altitude_km = 600;

% Pilot configuration
maxDelayspread_s = 3 * 30e-9;
maxDelayBins = ceil(maxDelayspread_s / delta_tau_s);
residualScatter_hz = 0.05 * fd_hz;
residualScatter_bins = ceil(residualScatter_hz / delta_nu_hz);

% QAM constellation (for reference / Fig 7)
qamAlphabet = qammod((0:ModOrder-1)', ModOrder, 'UnitAveragePower', true);

fprintf('=== LEO-NTN OTFS Full Paper Simulation ===\n');
fprintf('N=%d, M=%d, fd=%.0f Hz, delta_nu=%.1f Hz\n\n', ...
    numSC, ofdmSym, fd_hz, delta_nu_hz);

%% ========================================================================
%  Initialize LDPC Encoder/Decoder (DVB-S.2 standard)
%  ========================================================================
parityCheck_matrix = dvbs2ldpc(codeRate);
ldpcEncoder = comm.LDPCEncoder(parityCheck_matrix);
ldpcDecoder = comm.LDPCDecoder(parityCheck_matrix);
ldpcDecoder.MaximumIterationCount = maxLDPCIter;
noCodedbits = size(parityCheck_matrix, 2);          % Codeword length (64800)
infoBitsPerCW = noCodedbits * codeRate;              % Info bits per codeword (32400)

%% ========================================================================
%  Generate Common Data (single subframe)
%  ========================================================================
ddGridSize = [numSC, ofdmSym];
tfGridSize = zeros(numSC, ofdmSym);

% --- Uncoded data ---
txBits = randi([0, 1], numDC * ofdmSym * k, 1);
qamTx = qammod(txBits, ModOrder, 'InputType', 'bit', 'UnitAveragePower', true);
parallelTx = reshape(qamTx, [numDC, ofdmSym]);
guardbandTx = [zeros(1, ofdmSym); parallelTx];
guardbandTx = [guardbandTx(1:numDC/2, :); zeros(11, ofdmSym); ...
               guardbandTx(numDC/2+1:end, :)];

% --- Coded data (LDPC) ---
% Determine number of codewords that fit in one subframe
subframeBits = numDC * ofdmSym * k;
numCW = floor(subframeBits / noCodedbits);
if numCW < 1, numCW = 1; end
codedBitsTotal = numCW * noCodedbits;
padBits = max(0, subframeBits - codedBitsTotal);

% Generate info bits -> encode -> pad -> interleave -> QAM
infoBits_in = randi([0, 1], numCW * infoBitsPerCW, 1);
codedData_in = [];
for q = 1:numCW
    seg = infoBits_in((q-1)*infoBitsPerCW+1 : q*infoBitsPerCW);
    codedData_in = [codedData_in; ldpcEncoder(seg)];
end
paddedCoded = [codedData_in; zeros(padBits, 1)];
interleavedCoded = randintrlv(paddedCoded, 4831);
qamTx_coded = qammod(interleavedCoded, ModOrder, 'InputType', 'bit', 'UnitAveragePower', true);
parallelTx_coded = reshape(qamTx_coded, [numDC, ofdmSym]);
guardbandTx_coded = [zeros(1, ofdmSym); parallelTx_coded];
guardbandTx_coded = [guardbandTx_coded(1:numDC/2, :); zeros(11, ofdmSym); ...
                     guardbandTx_coded(numDC/2+1:end, :)];

% Pilot config
pilotCfg.kp = ceil(ofdmSym/2);
pilotCfg.lp = ceil(ddGridSize(1)/2);
pilotCfg.lGuard = 2 * maxDelayBins + 1;
pilotCfg.kGuard = 2 * residualScatter_bins + 2;
pilotCfg.maxDelayBins = maxDelayBins;
pilotCfg.maxDopplerBins = residualScatter_bins + 1;
pilotCfg.pilotBoostdB = 10;

%% ========================================================================
%  OTFS-Pilot Specific Data (DD-domain, independent of TF guard bands)
%  ========================================================================
numGuardTotal = (2*pilotCfg.kGuard+1) * (2*pilotCfg.lGuard+1);
numDataPosPilot = numSC * ofdmSym - numGuardTotal;

% Uncoded DD data
txBits_pilot = randi([0,1], numDataPosPilot * k, 1);
qamTx_pilot = qammod(txBits_pilot, ModOrder, 'InputType', 'bit', 'UnitAveragePower', true);

% Coded DD data (LDPC)
numCW_pilot = floor(numDataPosPilot * k / noCodedbits);
if numCW_pilot < 1, numCW_pilot = 1; end
codedBitsTotal_p = numCW_pilot * noCodedbits;
padBits_p = max(0, numDataPosPilot * k - codedBitsTotal_p);

infoBits_pilot = randi([0,1], numCW_pilot * infoBitsPerCW, 1);
codedData_p = [];
for q = 1:numCW_pilot
    seg = infoBits_pilot((q-1)*infoBitsPerCW+1 : q*infoBitsPerCW);
    codedData_p = [codedData_p; ldpcEncoder(seg)];
end
paddedCoded_p = [codedData_p; zeros(padBits_p, 1)];
interleavedCoded_p = randintrlv(paddedCoded_p, 4831);
qamTx_pilot_coded = qammod(interleavedCoded_p, ModOrder, 'InputType', 'bit', 'UnitAveragePower', true);

fprintf('  OTFS-Pilot: %d data positions, %d uncoded bits, %d coded CWs\n', ...
    numDataPosPilot, length(txBits_pilot), numCW_pilot);

%% ========================================================================
%  Fig 1: Three-Scheme BER vs Eb/No  (Monte Carlo averaged)
%  ========================================================================
fprintf('=== Fig 1: Three-Scheme BER vs Eb/No ===\n');
tFig = tic;

EbNo_fig1 = (0:2:24)';
numEbNo1 = length(EbNo_fig1);
numTrials1 = 50;   % Monte Carlo trials for statistical reliability

% SNR formulas (match main.m)
% OFDM/OTFS: guard-band overhead penalty
snr_ofdm = EbNo_fig1 + 10*log10(codeRate*k) + 10*log10(numDC/numSC);
snr_otfs = snr_ofdm;
% OTFS-Pilot: additional pilot+guard overhead penalty
snr_pilot = EbNo_fig1 + 10*log10(codeRate*k) + 10*log10(numDataPosPilot/(numSC*ofdmSym));

% Uncoded BER
ber1_ofdm = zeros(numEbNo1, 1);
ber1_otfs = zeros(numEbNo1, 1);
ber1_pilot = zeros(numEbNo1, 1);
% Coded BER (LDPC)
ber1_c_ofdm = zeros(numEbNo1, 1);
ber1_c_otfs = zeros(numEbNo1, 1);
ber1_c_pilot = zeros(numEbNo1, 1);

for m = 1:numEbNo1
    ber_o = 0; ber_t = 0; ber_p = 0;
    ber_co = 0; ber_ct = 0; ber_cp = 0;
    snr_o_m = snr_ofdm(m);  snr_t_m = snr_otfs(m);  snr_p_m = snr_pilot(m);
    nBitsTx = length(txBits);  nBitsInfo = length(infoBits_in);
    nBitsPilot = length(txBits_pilot);  nBitsInfoP = length(infoBits_pilot);

    parfor trial = 1:numTrials1
        % Local LDPC decoders (System Objects cannot cross worker boundary)
        locDec = comm.LDPCDecoder(parityCheck_matrix);
        locDec.MaximumIterationCount = maxLDPCIter;
        locDecP = comm.LDPCDecoder(parityCheck_matrix);
        locDecP.MaximumIterationCount = maxLDPCIter;

        % Generate new channel realization per trial
        [~, chI] = multipathChannel(cpSize, scs, tfGridSize, velocity, ntnConfig);
        % OTFS equalization channel: integer DD bins (matches applyChannelDD)
        chTF_eq = buildTF_AFC(chI, scs, cpSize, numSC, ofdmSym, fd_hz);
        % OFDM equalization channel: continuous delay/Doppler with sinc ICI loss
        chTF_ofdm = buildTF_OFDMeq(chI, scs, cpSize, numSC, ofdmSym, fd_hz);

        % ===================== OFDM =====================
        % Time-domain channel with per-sample Doppler (ICI)
        % Uncoded
        fTF = applyChannelTF_ICI(guardbandTx, chI, scs, cpSize, fd_hz);
        sp = mean(abs(fTF(:)).^2);
        nV = sp / 10^(snr_o_m/10);
        ns = sqrt(nV/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
        eq = conj(chTF_ofdm) ./ (abs(chTF_ofdm).^2 + nV) .* (fTF + ns);
        pR = eq; pR(numDC/2+1:numDC/2+11,:)=[]; pR(1,:)=[];
        rb = qamdemod(pR(:), ModOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
        ber_o = ber_o + sum(txBits ~= rb(1:nBitsTx)) / nBitsTx;

        % Coded
        fTF_c = applyChannelTF_ICI(guardbandTx_coded, chI, scs, cpSize, fd_hz);
        sp_c = mean(abs(fTF_c(:)).^2);
        nV_c = sp_c / 10^(snr_o_m/10);
        ns = sqrt(nV_c/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
        eq = conj(chTF_ofdm) ./ (abs(chTF_ofdm).^2 + nV_c) .* (fTF_c + ns);
        pR = eq; pR(numDC/2+1:numDC/2+11,:)=[]; pR(1,:)=[];
        % Empirical noise variance: captures ICI residual after one-tap equalization
        hardDec_o = qammod(qamdemod(pR(:), ModOrder, 'UnitAveragePower', true), ModOrder, 'UnitAveragePower', true);
        nV_c_emp = max(mean(abs(pR(:) - hardDec_o).^2), 1e-10);
        llr = qamdemod(pR(:), ModOrder, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', nV_c_emp);
        deintLLR = randdeintrlv(llr, 4831);
        deintLLR = deintLLR(1:numCW*noCodedbits);
        deintLLR = max(min(deintLLR, 50), -50);
        decBits = decodeLDPC_helper(deintLLR, locDec, numCW, noCodedbits);
        ber_co = ber_co + sum(infoBits_in ~= decBits(1:nBitsInfo)) / nBitsInfo;

        % ===================== OTFS =====================
        % Uncoded
        rxDD_t = applyChannelDD(guardbandTx, chI, scs, cpSize, fd_hz);
        fTF2 = ISFFT(rxDD_t);
        sp2 = mean(abs(fTF2(:)).^2);
        nV2 = sp2 / 10^(snr_t_m/10);
        ns = sqrt(nV2/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
        eq = conj(chTF_eq) ./ (abs(chTF_eq).^2 + nV2) .* (fTF2 + ns);
        rDD = SFFT(eq);
        pR = rDD; pR(numDC/2+1:numDC/2+11,:)=[]; pR(1,:)=[];
        rb = qamdemod(pR(:), ModOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
        ber_t = ber_t + sum(txBits ~= rb(1:nBitsTx)) / nBitsTx;

        % Coded
        rxDD_tc = applyChannelDD(guardbandTx_coded, chI, scs, cpSize, fd_hz);
        fTF2_c = ISFFT(rxDD_tc);
        sp2_c = mean(abs(fTF2_c(:)).^2);
        nV2_c = sp2_c / 10^(snr_t_m/10);
        ns = sqrt(nV2_c/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
        eq = conj(chTF_eq) ./ (abs(chTF_eq).^2 + nV2_c) .* (fTF2_c + ns);
        rDD = SFFT(eq);
        pR = rDD; pR(numDC/2+1:numDC/2+11,:)=[]; pR(1,:)=[];
        % Empirical noise variance: captures both thermal noise and MMSE-gain ISI
        hardDec = qammod(qamdemod(pR(:), ModOrder, 'UnitAveragePower', true), ModOrder, 'UnitAveragePower', true);
        nV2_emp = max(mean(abs(pR(:) - hardDec).^2), 1e-10);
        llr = qamdemod(pR(:), ModOrder, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', nV2_emp);
        deintLLR = randdeintrlv(llr, 4831);
        deintLLR = deintLLR(1:numCW*noCodedbits);
        deintLLR = max(min(deintLLR, 50), -50);
        decBits = decodeLDPC_helper(deintLLR, locDec, numCW, noCodedbits);
        ber_ct = ber_ct + sum(infoBits_in ~= decBits(1:nBitsInfo)) / nBitsInfo;

        % ===================== OTFS-Pilot =====================
        % Uncoded
        [dG, dI, ~, ~, pI] = pilotPatternDD(qamTx_pilot, ddGridSize(1), ofdmSym, pilotCfg);
        pI.maxDelayBins = pilotCfg.maxDelayBins;
        pI.maxDopplerBins = pilotCfg.maxDopplerBins;
        rF = applyChannelDD(dG, chI, scs, cpSize, fd_hz);
        sp3 = mean(abs(rF(:)).^2);
        nV3 = sp3 / 10^(snr_p_m/10);
        ns = sqrt(nV3/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
        rF_noisy = rF + ns;

        [~,~,~,nI] = ddChannelEstimate(rF_noisy, pI);
        pA = dG(pilotCfg.lp, pilotCfg.kp);
        pResp = zeros(ddGridSize);
        for pp = 1:nI.numPathsDetected
            lr = mod(pI.lp + nI.pathDelays(pp) - 1, ddGridSize(1)) + 1;
            kc = mod(pI.kp + nI.pathDopplers(pp) - 1, ddGridSize(2)) + 1;
            pResp(lr, kc) = nI.pathGains(pp);
        end
        cG = nI.pathGains / pA;
        estCh = struct('pathDelays_s', nI.pathDelays * delta_tau_s, ...
            'pathDopplers_Hz', nI.pathDopplers * delta_nu_hz + fd_hz, ...
            'pathGains', cG, 'numPaths', nI.numPathsDetected);
        chTF_est = buildTF_AFC(estCh, scs, cpSize, numSC, ofdmSym, fd_hz);
        rClean = rF_noisy - pResp;
        rClean_TF = ISFFT(rClean);
        eq_p = conj(chTF_est) ./ (abs(chTF_est).^2 + nV3) .* rClean_TF;
        rDD_eq = SFFT(eq_p);
        eqS = rDD_eq(dI);
        rb = qamdemod(eqS, ModOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
        ber_p = ber_p + sum(txBits_pilot ~= rb(1:nBitsPilot)) / nBitsPilot;

        % Coded
        [dG_c, dI_c, ~, ~, pI_c] = pilotPatternDD(qamTx_pilot_coded, ddGridSize(1), ofdmSym, pilotCfg);
        pI_c.maxDelayBins = pilotCfg.maxDelayBins;
        pI_c.maxDopplerBins = pilotCfg.maxDopplerBins;
        rF_c = applyChannelDD(dG_c, chI, scs, cpSize, fd_hz);
        sp3_c = mean(abs(rF_c(:)).^2);
        nV3_c = sp3_c / 10^(snr_p_m/10);
        ns = sqrt(nV3_c/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
        rF_c_noisy = rF_c + ns;

        [~,~,~,nI_c] = ddChannelEstimate(rF_c_noisy, pI_c);
        pA_c = dG_c(pilotCfg.lp, pilotCfg.kp);
        pResp_c = zeros(ddGridSize);
        for pp = 1:nI_c.numPathsDetected
            lr = mod(pI_c.lp + nI_c.pathDelays(pp) - 1, ddGridSize(1)) + 1;
            kc = mod(pI_c.kp + nI_c.pathDopplers(pp) - 1, ddGridSize(2)) + 1;
            pResp_c(lr, kc) = nI_c.pathGains(pp);
        end
        cG_c = nI_c.pathGains / pA_c;
        estCh_c = struct('pathDelays_s', nI_c.pathDelays * delta_tau_s, ...
            'pathDopplers_Hz', nI_c.pathDopplers * delta_nu_hz + fd_hz, ...
            'pathGains', cG_c, 'numPaths', nI_c.numPathsDetected);
        chTF_est_c = buildTF_AFC(estCh_c, scs, cpSize, numSC, ofdmSym, fd_hz);
        rClean_c = rF_c_noisy - pResp_c;
        rClean_c_TF = ISFFT(rClean_c);
        eq_pc = conj(chTF_est_c) ./ (abs(chTF_est_c).^2 + nV3_c) .* rClean_c_TF;
        rDD_eq_c = SFFT(eq_pc);
        eqS_c = rDD_eq_c(dI_c);
        % Empirical noise variance for OTFS-Pilot coded
        hardDec_p = qammod(qamdemod(eqS_c, ModOrder, 'UnitAveragePower', true), ModOrder, 'UnitAveragePower', true);
        nV3_c_emp = max(mean(abs(eqS_c - hardDec_p).^2), 1e-10);
        llr_p = qamdemod(eqS_c, ModOrder, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', nV3_c_emp);
        deintLLR_p = randdeintrlv(llr_p, 4831);
        deintLLR_p = deintLLR_p(1:numCW_pilot*noCodedbits);
        deintLLR_p = max(min(deintLLR_p, 50), -50);
        decBits_p = decodeLDPC_helper(deintLLR_p, locDecP, numCW_pilot, noCodedbits);
        ber_cp = ber_cp + sum(infoBits_pilot ~= decBits_p(1:nBitsInfoP)) / nBitsInfoP;
    end

    ber1_ofdm(m) = ber_o / numTrials1;
    ber1_otfs(m) = ber_t / numTrials1;
    ber1_pilot(m) = ber_p / numTrials1;
    ber1_c_ofdm(m) = ber_co / numTrials1;
    ber1_c_otfs(m) = ber_ct / numTrials1;
    ber1_c_pilot(m) = ber_cp / numTrials1;

    fprintf('  EbNo=%+3d dB: OFDM=%.4f/%.4f  OTFS=%.4f/%.4f  Pilot=%.4f/%.4f\n', ...
        EbNo_fig1(m), ber1_ofdm(m), ber1_c_ofdm(m), ...
        ber1_otfs(m), ber1_c_otfs(m), ber1_pilot(m), ber1_c_pilot(m));
end

% Plot Fig 1: 6 curves (3 uncoded solid + 3 coded dashed)
figure('Name', 'Fig1: Three-Scheme BER', 'Position', [100 100 750 580]);
% Uncoded (solid) — use offset markers to separate OFDM/OTFS visually
semilogy(EbNo_fig1, max(ber1_ofdm, 1e-6), 'g-o', 'LineWidth', 2, 'MarkerSize', 7); hold on;
semilogy(EbNo_fig1, max(ber1_otfs, 1e-6), 'r-s', 'LineWidth', 2, 'MarkerSize', 7);
semilogy(EbNo_fig1, max(ber1_pilot, 1e-6), 'b-d', 'LineWidth', 2, 'MarkerSize', 7);
% Coded (dashed) — use distinct marker sizes to separate C-OFDM/C-OTFS
semilogy(EbNo_fig1, max(ber1_c_ofdm, 1e-6), 'g--o', 'LineWidth', 1.5, 'MarkerSize', 7);
semilogy(EbNo_fig1, max(ber1_c_otfs, 1e-6), 'r--s', 'LineWidth', 1.5, 'MarkerSize', 5);
semilogy(EbNo_fig1, max(ber1_c_pilot, 1e-6), 'b--d', 'LineWidth', 1.5, 'MarkerSize', 7);
xlabel('E_b/N_0 (dB)'); ylabel('BER');
title(sprintf('Three-Scheme BER: %s, %s, Elev=%d° (%d trials)', ...
    ntnConfig.profile, ntnConfig.scenario, ntnConfig.elevAngle, numTrials1));
legend('OFDM (TF-MMSE)', 'OTFS (TF-MMSE)', 'OTFS-Pilot (Est+TF-MMSE)', ...
    'C-OFDM (LDPC+TF-MMSE)', 'C-OTFS (LDPC+TF-MMSE)', 'C-OTFS-Pilot (LDPC+Est+TF-MMSE)', ...
    'Location', 'southwest', 'FontSize', 9);
grid on; ylim([1e-6, 1]);
set(gca, 'FontSize', 11);

% No annotations — Monte Carlo averaging reveals OTFS diversity advantage

saveas(gcf, 'fig1_three_scheme_ber.fig');
print(gcf, 'fig1_three_scheme_ber', '-dpng', '-r150');
fprintf('  Fig 1 done (%.1f s)\n\n', toc(tFig));


%% ========================================================================
%  Fig 4: Pilot Overhead vs Pre-compensation Ratio (Analytical)
%  ========================================================================
fprintf('=== Fig 4: Pilot Overhead ===\n');
tFig = tic;

precompRatios = 0:0.05:1.0;
numRatios = length(precompRatios);
overheadPct = zeros(numRatios, 1);
kGuardVals = zeros(numRatios, 1);
lGuardFixed = 2 * maxDelayBins + 1;

for idx = 1:numRatios
    ratio = precompRatios(idx);
    res_hz = (1 - ratio) * fd_hz + 0.05 * fd_hz;
    res_bins = ceil(res_hz / delta_nu_hz);
    kG = min(2 * res_bins + 2, floor(ofdmSym/2) - 1);
    kGuardVals(idx) = kG;
    overheadPct(idx) = 100 * (2*kG+1) * (2*lGuardFixed+1) / (numSC * ofdmSym);
end

feasMask = kGuardVals < (floor(ofdmSym/2) - 1);
dataCap = 1 - overheadPct / 100;

% Find feasibility threshold
feasThresh = NaN;
for idx = 1:numRatios
    if feasMask(idx), feasThresh = precompRatios(idx)*100; break; end
end

figure('Name', 'Fig4: Pilot Overhead', 'Position', [100 100 900 400]);

subplot(1,2,1);
plot(precompRatios*100, overheadPct, 'b-o', 'LineWidth', 2, 'MarkerSize', 4);
hold on;
if ~isnan(feasThresh)
    xline(feasThresh, 'r--', sprintf('Feasibility (%.0f%%)', feasThresh), ...
        'LineWidth', 1.5, 'LabelOrientation', 'aligned', 'LabelVerticalAlignment', 'bottom');
    patch([0 feasThresh feasThresh 0], [0 0 100 100], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end
xlabel('Doppler Pre-comp Ratio (%)'); ylabel('Overhead (%)');
title('(a) Pilot + Guard Overhead');
grid on; ylim([0 100]); set(gca, 'FontSize', 11);

subplot(1,2,2);
yyaxis left;
barIdx = [];
for t = [0 20 50 70 80 90 100]
    [~,bi] = min(abs(precompRatios*100-t)); barIdx = [barIdx, bi];
end
bar(precompRatios(barIdx)*100, kGuardVals(barIdx), 0.5);
ylabel('k_{Guard} (bins)');
yyaxis right;
plot(precompRatios*100, dataCap*100, 'r-s', 'LineWidth', 2, 'MarkerSize', 4);
ylabel('Data Capacity (%)');
xlabel('Pre-comp Ratio (%)');
title('(b) Guard Size & Capacity');
legend('k_{Guard}', 'Capacity', 'Location', 'east');
grid on; set(gca, 'FontSize', 11);

saveas(gcf, 'fig4_overhead_vs_precomp.fig');
print(gcf, 'fig4_overhead_vs_precomp', '-dpng', '-r150');
fprintf('  Fig 4 done (%.1f s)\n\n', toc(tFig));

%% ========================================================================
%  Fig 5: Velocity & Range RMSE vs SNR — OTFS-Pilot vs OFDM-Pilot vs CRB
%  ========================================================================
fprintf('=== Fig 5: Navigation Accuracy vs SNR ===\n');
tFig = tic;

EbNo_nav = (0:2:24)';
numNavEbNo = length(EbNo_nav);
numNavTrials = 500;
snr_nav = EbNo_nav + 10*log10(codeRate*k) + 10*log10(numDC/numSC) ;

velRMSE = zeros(numNavEbNo, 1);
rangeRMSE = zeros(numNavEbNo, 1);
velRMSE_ofdm = zeros(numNavEbNo, 1);
rangeRMSE_ofdm = zeros(numNavEbNo, 1);

for m = 1:numNavEbNo
    vErr2 = zeros(numNavTrials, 1);
    rErr2 = zeros(numNavTrials, 1);
    vErr2_o = zeros(numNavTrials, 1);
    rErr2_o = zeros(numNavTrials, 1);

    snr_nav_m = snr_nav(m);
    parfor trial = 1:numNavTrials
        [~, chN] = multipathChannel(cpSize, scs, tfGridSize, velocity, ntnConfig);

        % ---- Add random residual propagation delay ----
        bulkDelayBins = 2 + 6*rand();
        bulkDelay_s   = bulkDelayBins * delta_tau_s;
        chN.pathDelays_s = chN.pathDelays_s + bulkDelay_s;

        % True LoS (always first tap by construction in multipathChannel)
        trueVel = chN.pathDopplers_Hz(1) * c_light / fc;
        trueRange = chN.pathDelays_s(1) * c_light;

        % ============ OTFS-Pilot sensing ============
        [dG, ~, ~, ~, pI] = pilotPatternDD(qamTx_pilot, ddGridSize(1), ofdmSym, pilotCfg);
        pI.maxDelayBins = 10;   % Cover bulk delay (2-8 bins) + multipath
        pI.maxDopplerBins = pilotCfg.maxDopplerBins;
        rF = applyChannelDD(dG, chN, scs, cpSize, fd_hz, true);

        sp = mean(abs(rF(:)).^2);
        nV = sp * 10^(-snr_nav_m/10);
        rF = rF + sqrt(nV/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));

        [~,~,~,nE] = ddChannelEstimate(rF, pI);

        estDoppler_Hz = nE.losFracDoppler * delta_nu_hz + fd_hz;
        estDelay_s = nE.losFracDelay * delta_tau_s;
        estVel = estDoppler_Hz * c_light / fc;
        estRange = estDelay_s * c_light;

        vErr2(trial) = (estVel - trueVel)^2;
        rErr2(trial) = (estRange - trueRange)^2;

        % ============ OFDM-Pilot sensing (LS estimation) ============
        % Realistic OFDM: send known signal through TF-ICI channel, LS estimate
        fTF_nav = applyChannelTF_ICI(guardbandTx, chN, scs, cpSize, fd_hz);
        sp_o = mean(abs(fTF_nav(:)).^2);
        nV_o = sp_o * 10^(-snr_nav_m/10);
        fTF_nav_noisy = fTF_nav + sqrt(nV_o/2) * ...
            (randn(ddGridSize) + 1j*randn(ddGridSize));

        % LS channel estimation: H_est = Y/X at data positions only
        H_est = zeros(ddGridSize);
        dataMask_nav = abs(guardbandTx) > 1e-10;
        H_est(dataMask_nav) = fTF_nav_noisy(dataMask_nav) ./ guardbandTx(dataMask_nav);

        % Interpolate guard band positions (DC + center guard) per OFDM symbol
        guardRows = find(~dataMask_nav(:,1));   % Same pattern for all symbols
        dataRows  = find(dataMask_nav(:,1));
        for sym = 1:ofdmSym
            H_est(guardRows, sym) = interp1(dataRows, H_est(dataRows, sym), ...
                guardRows, 'linear', 'extrap');
        end

        % Convert estimated TF channel to DD domain
        % Note: applyChannelTF_ICI already performs CFO correction internally,
        % so H_est already has bulk Doppler removed. No additional pre-comp needed.
        ddResp = fftshift(fft(ifft(H_est, [], 1), [], 2));

        % Search for LoS peak in RESTRICTED delay-Doppler region
        % (prevents guard-band spectral leakage artifacts from being selected)
        delayCentre = floor(numSC/2) + 1;
        dopplerCentre = floor(ofdmSym/2) + 1;
        ofdmScanRows = delayCentre : min(numSC, delayCentre + 20);
        ofdmScanCols = max(1, dopplerCentre - 10) : min(ofdmSym, dopplerCentre + 10);
        ddResp_scan = abs(ddResp(ofdmScanRows, ofdmScanCols));
        [~, scanPeakIdx] = max(ddResp_scan(:));
        [scanRow, scanCol] = ind2sub(size(ddResp_scan), scanPeakIdx);
        peakDelay = ofdmScanRows(scanRow);
        peakDoppler = ofdmScanCols(scanCol);

        % Convert peak indices to physical units
        estDelay_ofdm_s = (peakDelay - delayCentre) * delta_tau_s;
        estResidualDoppler_Hz = (peakDoppler - dopplerCentre) * delta_nu_hz;

        % Quinn fractional refinement — Doppler dimension
        if peakDoppler > 1 && peakDoppler < ofdmSym
            Xm1 = ddResp(peakDelay, peakDoppler-1);
            X0  = ddResp(peakDelay, peakDoppler);
            Xp1 = ddResp(peakDelay, peakDoppler+1);
            ap = real(Xp1/X0); am = real(Xm1/X0);
            dp = -ap/(1-ap); dm = am/(1-am);
            if abs(dp) < abs(dm), dk = dp; else, dk = dm; end
            dk = max(-0.5, min(0.5, dk));
            estResidualDoppler_Hz = (peakDoppler - dopplerCentre + dk) * delta_nu_hz;
        end
        % Quinn fractional refinement — Delay dimension
        if peakDelay > 1 && peakDelay < numSC
            Xm1 = ddResp(peakDelay-1, peakDoppler);
            X0  = ddResp(peakDelay, peakDoppler);
            Xp1 = ddResp(peakDelay+1, peakDoppler);
            ap = real(Xp1/X0); am = real(Xm1/X0);
            dp = -ap/(1-ap); dm = am/(1-am);
            if abs(dp) < abs(dm), dl = dp; else, dl = dm; end
            dl = max(-0.5, min(0.5, dl));
            estDelay_ofdm_s = (peakDelay - delayCentre + dl) * delta_tau_s;
        end
        % Add back pre-compensated bulk Doppler
        estDoppler_ofdm_Hz = estResidualDoppler_Hz + fd_hz;
        estVel_o = estDoppler_ofdm_Hz * c_light / fc;
        estRange_o = abs(estDelay_ofdm_s) * c_light;

        vErr2_o(trial) = (estVel_o - trueVel)^2;
        rErr2_o(trial) = (estRange_o - trueRange)^2;
    end

    % Median-based RMSE (robust to outlier trials from fading/detection failures)
    % median(X^2) is a robust scale estimator; sqrt gives "median RMSE"
    velRMSE(m) = sqrt(median(vErr2));
    rangeRMSE(m) = sqrt(median(rErr2));
    velRMSE_ofdm(m) = sqrt(median(vErr2_o));
    rangeRMSE_ofdm(m) = sqrt(median(rErr2_o));
    fprintf('  EbNo=%+3d: OTFS vel=%.2f/rng=%.3f  OFDM vel=%.2f/rng=%.3f\n', ...
        EbNo_nav(m), velRMSE(m), rangeRMSE(m), velRMSE_ofdm(m), rangeRMSE_ofdm(m));
end

% CRB — Kay's frequency estimation CRB, separate for each scheme
% σ_f = √12 / (2π·Ts·√(SNR_eff · M · (M²-1)))
% Velocity: σ_v = (c/fc) · σ_f
% Range:    σ_r = c · σ_τ, σ_τ = √12 / (2π·Δf·√(SNR_eff · N · (N²-1)))
snr_lin = 10.^(snr_nav/10);
pilotBoost_lin = 10^(pilotCfg.pilotBoostdB/10);

% OFDM CRB: data-aided sensing, coherent over N×M observations
% Velocity: N subcarriers boost Doppler SNR by factor N
crb_vel_ofdm = (c_light/fc) * sqrt(12) ./ ...
    (2*pi*Ts_sym * sqrt(numSC .* snr_lin * ofdmSym * (ofdmSym^2-1)));
% Range: M symbols boost delay SNR by factor M
crb_range_ofdm = c_light * sqrt(12) ./ ...
    (2*pi*scs * sqrt(ofdmSym .* snr_lin * numSC * (numSC^2-1)));

% OTFS CRB: single pilot with boost, M Doppler / N delay observations
crb_vel_otfs = (c_light/fc) * sqrt(12) ./ ...
    (2*pi*Ts_sym * sqrt(pilotBoost_lin .* snr_lin * ofdmSym * (ofdmSym^2-1)));
crb_range_otfs = c_light * sqrt(12) ./ ...
    (2*pi*scs * sqrt(pilotBoost_lin .* snr_lin * numSC * (numSC^2-1)));

figure('Name', 'Fig5: Nav Accuracy vs SNR', 'Position', [100 100 900 400]);

subplot(1,2,1);
semilogy(EbNo_nav, velRMSE_ofdm, 'g-o', 'LineWidth', 2, 'MarkerSize', 6); hold on;
semilogy(EbNo_nav, velRMSE, 'b-d', 'LineWidth', 2, 'MarkerSize', 7);
semilogy(EbNo_nav, crb_vel_ofdm, 'g--', 'LineWidth', 1.5);
semilogy(EbNo_nav, crb_vel_otfs, 'b--', 'LineWidth', 1.5);
xlabel('E_b/N_0 (dB)'); ylabel('Velocity RMSE (m/s)');
title('(a) Velocity Estimation');
legend('OFDM (data-aided)', 'OTFS DD Pilot', ...
    'CRB_{OFDM}', 'CRB_{OTFS}', 'Location', 'southwest');
grid on; set(gca, 'FontSize', 11);

subplot(1,2,2);
semilogy(EbNo_nav, rangeRMSE_ofdm, 'g-o', 'LineWidth', 2, 'MarkerSize', 6); hold on;
semilogy(EbNo_nav, rangeRMSE, 'b-d', 'LineWidth', 2, 'MarkerSize', 7);
semilogy(EbNo_nav, crb_range_ofdm, 'g--', 'LineWidth', 1.5);
semilogy(EbNo_nav, crb_range_otfs, 'b--', 'LineWidth', 1.5);
xlabel('E_b/N_0 (dB)'); ylabel('Range RMSE (m)');
title('(b) Range Estimation');
legend('OFDM (data-aided)', 'OTFS DD Pilot', ...
    'CRB_{OFDM}', 'CRB_{OTFS}', 'Location', 'southwest');
grid on; set(gca, 'FontSize', 11);

saveas(gcf, 'fig5_nav_accuracy_snr.fig');
print(gcf, 'fig5_nav_accuracy_snr', '-dpng', '-r150');
fprintf('  Fig 5 done (%.1f s)\n\n', toc(tFig));


%% ========================================================================
%  Fig 7: DD Grid Pilot Pattern Visualization
%  ========================================================================
fprintf('=== Fig 7: DD Grid Visualization ===\n');
tFig = tic;

figure('Name', 'Fig7: DD Grid Comparison', 'Position', [100 100 1000 450]);

% (a) WITHOUT pre-comp
fullDop_bins = ceil(fd_hz / delta_nu_hz);
kG_noPre = min(2*fullDop_bins+2, floor(ofdmSym/2)-1);

subplot(1,2,1);
if kG_noPre >= floor(ofdmSym/2) - 1
    imagesc(zeros(ddGridSize)); colormap(hot);
    text(ofdmSym/2, ddGridSize(1)/2, ...
        sprintf('INFEASIBLE\nk_{Guard}=%d > M/2=%d\nOverhead > 100%%', ...
        kG_noPre, floor(ofdmSym/2)), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, ...
        'Color', 'w', 'FontWeight', 'bold', 'BackgroundColor', [0.8 0 0]);
    title(sprintf('(a) No Pre-comp (k_{Guard}=%d)', kG_noPre));
else
    pcNP = pilotCfg;
    pcNP.kGuard = kG_noPre;
    [~,dI_np,pI_np,gI_np,info_np] = pilotPatternDD(qamTx_pilot, ddGridSize(1), ofdmSym, pcNP);
    vis = zeros(ddGridSize); vis(dI_np)=0.3; vis(gI_np)=0; vis(pI_np)=1;
    imagesc(vis); colormap(hot);
    title(sprintf('(a) No Pre-comp: Overhead=%.1f%%', info_np.overheadPercent));
    hold on;
    rectangle('Position', [pcNP.kp-pcNP.kGuard-0.5, pcNP.lp-pcNP.lGuard-0.5, ...
        2*pcNP.kGuard+1, 2*pcNP.lGuard+1], 'EdgeColor','c','LineWidth',2,'LineStyle','--');
    hold off;
end
xlabel('Doppler Index'); ylabel('Delay Index');
set(gca, 'FontSize', 11);

% (b) WITH pre-comp
[~,dI_wp,pI_wp,gI_wp,info_wp] = pilotPatternDD(qamTx_pilot, ddGridSize(1), ofdmSym, pilotCfg);
subplot(1,2,2);
vis = zeros(ddGridSize); vis(dI_wp)=0.3; vis(gI_wp)=0; vis(pI_wp)=1;
imagesc(vis); colormap(hot);
title(sprintf('(b) With Pre-comp: Overhead=%.1f%%', info_wp.overheadPercent));
xlabel('Doppler Index'); ylabel('Delay Index');
hold on;
rectangle('Position', [pilotCfg.kp-pilotCfg.kGuard-0.5, pilotCfg.lp-pilotCfg.lGuard-0.5, ...
    2*pilotCfg.kGuard+1, 2*pilotCfg.lGuard+1], 'EdgeColor','c','LineWidth',2,'LineStyle','--');
hold off;
set(gca, 'FontSize', 11);

saveas(gcf, 'fig7_dd_grid_comparison.fig');
print(gcf, 'fig7_dd_grid_comparison', '-dpng', '-r150');
fprintf('  Fig 7 done (%.1f s)\n\n', toc(tFig));

%% ========================================================================
%  Summary Table
%  ========================================================================
fprintf('\n============================================================\n');
fprintf('  PAPER SIMULATION SUMMARY\n');
fprintf('============================================================\n');
fprintf('  System: N=%d, M=%d, %d-QAM, SCS=%.0f kHz\n', numSC, ofdmSym, ModOrder, scs/1e3);
fprintf('  Channel: 3GPP TR 38.811 NTN-TDL-C\n');
fprintf('  LEO: alt=%d km, v=%d km/h, fc=%.1f GHz\n', ...
    ntnConfig.altitude_km, velocity, fc/1e9);
fprintf('  Max Doppler: %.0f Hz (%.0f bins)\n', fd_hz, fd_hz/delta_nu_hz);
fprintf('  Pilot overhead: %.1f%% (with pre-comp)\n', info_wp.overheadPercent);
fprintf('\n  --- Fig 1: BER at Eb/No=20 dB ---\n');
idx20 = find(EbNo_fig1==20, 1);
if ~isempty(idx20)
    fprintf('    OFDM:       %.2e\n', ber1_ofdm(idx20));
    fprintf('    OTFS:       %.2e\n', ber1_otfs(idx20));
    fprintf('    OTFS-Pilot: %.2e\n', ber1_pilot(idx20));
end
fprintf('\n  --- Fig 5: Nav at Eb/No=20 dB ---\n');
idx20n = find(EbNo_nav==20, 1);
if ~isempty(idx20n)
    fprintf('    Velocity RMSE: %.2f m/s (%.2f km/h)\n', velRMSE(idx20n), velRMSE(idx20n)*3.6);
    fprintf('    Range RMSE:    %.2f m\n', rangeRMSE(idx20n));
end
fprintf('============================================================\n');
fprintf('\n4 figures generated. Total time: %.1f s (%.1f min)\n', toc(tTotal), toc(tTotal)/60);
