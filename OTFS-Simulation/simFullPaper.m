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
%  Generates 7 Figures:
%    Fig 1: Three-scheme BER vs Eb/No (baseline)
%    Fig 2: BER vs Elevation Angle
%    Fig 3: BER vs Propagation Scenario
%    Fig 4: Pilot Overhead vs Pre-compensation Ratio (analytical)
%    Fig 5: Velocity/Range RMSE vs SNR (with CRB)
%    Fig 6: Navigation Accuracy vs Elevation Angle
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
%  Fig 2: BER vs Elevation Angle
%  ========================================================================
fprintf('=== Fig 2: BER vs Elevation Angle ===\n');
tFig = tic;

elevAngles = [10, 20, 30, 50, 70, 90];
numElev = length(elevAngles);
fixedEbNo = [12, 16];                  % Two Eb/No reference points
numFixedEbNo = length(fixedEbNo);
numTrials2 = 20;

ber2_ofdm  = zeros(numElev, numFixedEbNo);
ber2_otfs  = zeros(numElev, numFixedEbNo);
ber2_pilot = zeros(numElev, numFixedEbNo);

for ei = 1:numElev
    ntnCfg_e = ntnConfig;
    ntnCfg_e.elevAngle = elevAngles(ei);

    for fi = 1:numFixedEbNo
        ebno = fixedEbNo(fi);
        snr_o = ebno + 10*log10(codeRate*k) + 10*log10(numDC/numSC);
        snr_t = snr_o;
        snr_p = ebno + 10*log10(codeRate*k) + 10*log10(numDataPosPilot/(numSC*ofdmSym));

        nBitsTx2 = length(txBits);
        nBitsPilot2 = length(txBits_pilot);
        snr_o_fi = snr_o; snr_t_fi = snr_t; snr_p_fi = snr_p;
        ber_o_acc = 0; ber_t_acc = 0; ber_p_acc = 0;
        parfor trial = 1:numTrials2
            [~, chI] = multipathChannel(cpSize, scs, tfGridSize, velocity, ntnCfg_e);

            % Equalization channels
            chTF_eq = buildTF_AFC(chI, scs, cpSize, numSC, ofdmSym, fd_hz);
            chTF_ofdm = buildTF_OFDMeq(chI, scs, cpSize, numSC, ofdmSym, fd_hz);

            % OFDM (time-domain channel with ICI)
            fTF = applyChannelTF_ICI(guardbandTx, chI, scs, cpSize, fd_hz);
            sp = mean(abs(fTF(:)).^2);
            nV = sp / 10^(snr_o_fi/10);
            ns = sqrt(nV/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
            eq = conj(chTF_ofdm) ./ (abs(chTF_ofdm).^2 + nV) .* (fTF + ns);
            pR = eq; pR(numDC/2+1:numDC/2+11,:)=[]; pR(1,:)=[];
            rb = qamdemod(pR(:), ModOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
            ber_o_acc = ber_o_acc + sum(txBits ~= rb(1:nBitsTx2)) / nBitsTx2;

            % OTFS (DD -> channel -> DD -> ISFFT -> TF-MMSE -> SFFT)
            rxDD_t = applyChannelDD(guardbandTx, chI, scs, cpSize, fd_hz);
            fTF2 = ISFFT(rxDD_t);
            sp2 = mean(abs(fTF2(:)).^2);
            nV2 = sp2 / 10^(snr_t_fi/10);
            ns2 = sqrt(nV2/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
            eq2 = conj(chTF_eq) ./ (abs(chTF_eq).^2 + nV2) .* (fTF2 + ns2);
            rDD = SFFT(eq2);
            pR = rDD; pR(numDC/2+1:numDC/2+11,:)=[]; pR(1,:)=[];
            rb = qamdemod(pR(:), ModOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
            ber_t_acc = ber_t_acc + sum(txBits ~= rb(1:nBitsTx2)) / nBitsTx2;

            % OTFS-Pilot
            [dG, dI, ~, ~, pI] = pilotPatternDD(qamTx_pilot, ddGridSize(1), ofdmSym, pilotCfg);
            pI.maxDelayBins = pilotCfg.maxDelayBins;
            pI.maxDopplerBins = pilotCfg.maxDopplerBins;
            rF = applyChannelDD(dG, chI, scs, cpSize, fd_hz);
            sp3 = mean(abs(rF(:)).^2);
            nV3 = sp3 / 10^(snr_p_fi/10);
            ns3 = sqrt(nV3/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
            rF_noisy = rF + ns3;

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
            ber_p_acc = ber_p_acc + sum(txBits_pilot ~= rb(1:nBitsPilot2)) / nBitsPilot2;
        end
        ber2_ofdm(ei, fi)  = ber_o_acc / numTrials2;
        ber2_otfs(ei, fi)  = ber_t_acc / numTrials2;
        ber2_pilot(ei, fi) = ber_p_acc / numTrials2;

        fprintf('  Elev=%2d° EbNo=%2d: OFDM=%.4f OTFS=%.4f Pilot=%.4f\n', ...
            elevAngles(ei), ebno, ber2_ofdm(ei,fi), ber2_otfs(ei,fi), ber2_pilot(ei,fi));
    end
end

% Plot Fig 2
figure('Name', 'Fig2: BER vs Elevation', 'Position', [100 100 900 400]);
for fi = 1:numFixedEbNo
    subplot(1, numFixedEbNo, fi);
    semilogy(elevAngles, max(ber2_ofdm(:,fi), 1e-6), 'g-o', 'LineWidth', 2, 'MarkerSize', 7);
    hold on;
    semilogy(elevAngles, max(ber2_otfs(:,fi), 1e-6), 'r-s', 'LineWidth', 2, 'MarkerSize', 7);
    semilogy(elevAngles, max(ber2_pilot(:,fi), 1e-6), 'b-d', 'LineWidth', 2, 'MarkerSize', 7);
    xlabel('Elevation Angle (°)'); ylabel('BER');
    title(sprintf('E_b/N_0 = %d dB', fixedEbNo(fi)));
    legend('OFDM', 'OTFS', 'OTFS-Pilot', 'Location', 'northeast');
    grid on; ylim([1e-5, 1]);
    set(gca, 'FontSize', 11);
    xticks(elevAngles);
end
sgtitle('BER vs Satellite Elevation Angle (NTN-TDL-C, Suburban)');
saveas(gcf, 'fig2_ber_vs_elevation.fig');
print(gcf, 'fig2_ber_vs_elevation', '-dpng', '-r150');
fprintf('  Fig 2 done (%.1f s)\n\n', toc(tFig));

%% ========================================================================
%  Fig 3: BER vs Propagation Scenario
%  ========================================================================
fprintf('=== Fig 3: BER vs Scenario ===\n');
tFig = tic;

scenarios = {'DenseUrban', 'Urban', 'Suburban', 'Rural'};
numScen = length(scenarios);
numTrials3 = 20;
fixedEbNo3 = 15;   % Single reference point

snr_o3 = fixedEbNo3 + 10*log10(codeRate*k) + 10*log10(numDC/numSC);
snr_t3 = snr_o3;
snr_p3 = fixedEbNo3 + 10*log10(codeRate*k) + 10*log10(numDataPosPilot/(numSC*ofdmSym));

ber3_ofdm  = zeros(numScen, 1);
ber3_otfs  = zeros(numScen, 1);
ber3_pilot = zeros(numScen, 1);

for si = 1:numScen
    ntnCfg_s = ntnConfig;
    ntnCfg_s.scenario = scenarios{si};

    nBitsTx3 = length(txBits);
    nBitsPilot3 = length(txBits_pilot);
    ber_o_acc = 0; ber_t_acc = 0; ber_p_acc = 0;
    parfor trial = 1:numTrials3
        [~, chI] = multipathChannel(cpSize, scs, tfGridSize, velocity, ntnCfg_s);

        % Equalization channels
        chTF_eq = buildTF_AFC(chI, scs, cpSize, numSC, ofdmSym, fd_hz);
        chTF_ofdm = buildTF_OFDMeq(chI, scs, cpSize, numSC, ofdmSym, fd_hz);

        % OFDM (time-domain channel with ICI)
        fTF = applyChannelTF_ICI(guardbandTx, chI, scs, cpSize, fd_hz);
        sp = mean(abs(fTF(:)).^2);
        nV = sp / 10^(snr_o3/10);
        ns = sqrt(nV/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
        eq = conj(chTF_ofdm) ./ (abs(chTF_ofdm).^2 + nV) .* (fTF + ns);
        pR = eq; pR(numDC/2+1:numDC/2+11,:)=[]; pR(1,:)=[];
        rb = qamdemod(pR(:), ModOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
        ber_o_acc = ber_o_acc + sum(txBits ~= rb(1:nBitsTx3)) / nBitsTx3;

        % OTFS (DD -> channel -> DD -> ISFFT -> TF-MMSE -> SFFT)
        rxDD_t = applyChannelDD(guardbandTx, chI, scs, cpSize, fd_hz);
        fTF2 = ISFFT(rxDD_t);
        sp2 = mean(abs(fTF2(:)).^2);
        nV2 = sp2 / 10^(snr_t3/10);
        ns2 = sqrt(nV2/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
        eq2 = conj(chTF_eq) ./ (abs(chTF_eq).^2 + nV2) .* (fTF2 + ns2);
        rDD = SFFT(eq2);
        pR = rDD; pR(numDC/2+1:numDC/2+11,:)=[]; pR(1,:)=[];
        rb = qamdemod(pR(:), ModOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
        ber_t_acc = ber_t_acc + sum(txBits ~= rb(1:nBitsTx3)) / nBitsTx3;

        % OTFS-Pilot
        [dG, dI, ~, ~, pI] = pilotPatternDD(qamTx_pilot, ddGridSize(1), ofdmSym, pilotCfg);
        pI.maxDelayBins = pilotCfg.maxDelayBins;
        pI.maxDopplerBins = pilotCfg.maxDopplerBins;
        rF = applyChannelDD(dG, chI, scs, cpSize, fd_hz);
        sp3 = mean(abs(rF(:)).^2);
        nV3 = sp3 / 10^(snr_p3/10);
        ns3 = sqrt(nV3/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
        rF_noisy = rF + ns3;

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
        ber_p_acc = ber_p_acc + sum(txBits_pilot ~= rb(1:nBitsPilot3)) / nBitsPilot3;
    end
    ber3_ofdm(si)  = ber_o_acc / numTrials3;
    ber3_otfs(si)  = ber_t_acc / numTrials3;
    ber3_pilot(si) = ber_p_acc / numTrials3;

    fprintf('  [%s] OFDM=%.4f  OTFS=%.4f  Pilot=%.4f\n', ...
        scenarios{si}, ber3_ofdm(si), ber3_otfs(si), ber3_pilot(si));
end

% Plot Fig 3 — Grouped bar chart with data labels
figure('Name', 'Fig3: BER vs Scenario', 'Position', [100 100 780 500]);

barData = [ber3_ofdm, ber3_otfs, ber3_pilot];
x3 = 1:numScen;
b = bar(x3, barData, 0.82, 'grouped');

% Modern color scheme
colors = [0.35 0.70 0.35;    % OFDM  - green
          0.85 0.33 0.10;    % OTFS  - orange-red
          0.20 0.40 0.85];   % Pilot - blue
for bi = 1:3
    b(bi).FaceColor = 'flat';
    b(bi).CData = repmat(colors(bi,:), numScen, 1);
    b(bi).EdgeColor = colors(bi,:) * 0.6;
    b(bi).LineWidth = 1.0;
end

set(gca, 'YScale', 'log');
set(gca, 'XTick', x3, 'XTickLabel', ...
    {'Dense Urban', 'Urban', 'Suburban', 'Rural'});
set(gca, 'FontSize', 11, 'FontName', 'Arial');
set(gca, 'TickDir', 'out', 'Box', 'off');

% Y-axis: tighter range around data
allBER = barData(:);
yLo = 10^(floor(log10(min(allBER(allBER>0)))-0.5));
yHi = 10^(ceil(log10(max(allBER))+0.3));
ylim([yLo, yHi]);

ylabel('BER'); xlabel('Propagation Scenario');
title(sprintf('BER by Scenario (E_b/N_0 = %d dB, Elev = %d\\circ)', ...
    fixedEbNo3, ntnConfig.elevAngle), 'FontSize', 13);

% Add data labels on top of each bar
for bi = 1:3
    xData = b(bi).XEndPoints;
    yData = b(bi).YEndPoints;
    for j = 1:numScen
        if yData(j) > 0
            expVal = floor(log10(yData(j)));
            mantissa = yData(j) / 10^expVal;
            if abs(mantissa - 1) < 0.05
                labelStr = sprintf('10^{%d}', expVal);
            else
                labelStr = sprintf('%.1f\\times10^{%d}', mantissa, expVal);
            end
            text(xData(j), yData(j)*1.4, labelStr, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontSize', 8, 'Color', colors(bi,:)*0.7, 'FontWeight', 'bold');
        end
    end
end

legend('OFDM (TF-MMSE)', 'OTFS (TF-MMSE)', 'OTFS-Pilot (Est+TF-MMSE)', ...
    'Location', 'northeast', 'FontSize', 10);
grid on;
set(gca, 'GridAlpha', 0.15, 'MinorGridLineStyle', 'none');

saveas(gcf, 'fig3_ber_vs_scenario.fig');
print(gcf, 'fig3_ber_vs_scenario', '-dpng', '-r150');
fprintf('  Fig 3 done (%.1f s)\n\n', toc(tFig));

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
numNavTrials = 200;
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
        [chTF_nav, chN] = multipathChannel(cpSize, scs, tfGridSize, velocity, ntnConfig);

        % ---- Add random residual propagation delay ----
        bulkDelayBins = 2 + 6*rand();
        bulkDelay_s   = bulkDelayBins * delta_tau_s;
        chN.pathDelays_s = chN.pathDelays_s + bulkDelay_s;

        % Apply corresponding phase ramp to TF channel matrix
        sc_idx = (0:numSC-1)';
        phaseRamp = exp(-2j*pi * sc_idx * scs * bulkDelay_s);
        chTF_nav = chTF_nav .* phaseRamp;

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

        % ============ OFDM-Pilot sensing ============
        sp_o = mean(abs(chTF_nav(:)).^2);
        nV_o = sp_o * 10^(-snr_nav_m/10);
        chTF_noisy = chTF_nav + sqrt(nV_o/2) * ...
            (randn(size(chTF_nav)) + 1j*randn(size(chTF_nav)));

        % Doppler pre-compensation in TF domain: remove bulk fd
        m_precomp = exp(-1j*2*pi*fd_hz*(1:ofdmSym)*Ts_sym);
        chTF_precomp = chTF_noisy .* m_precomp;

        % Symplectic FFT: TF → delay-Doppler response
        ddResp = fftshift(fft(ifft(chTF_precomp, [], 1), [], 2));

        % Search for LoS peak in delay-Doppler map
        [~, peakIdx] = max(abs(ddResp(:)));
        [peakDelay, peakDoppler] = ind2sub(size(ddResp), peakIdx);

        % Convert peak indices to physical units
        delayCentre = floor(numSC/2) + 1;
        estDelay_ofdm_s = (peakDelay - delayCentre) * delta_tau_s;

        dopplerCentre = floor(ofdmSym/2) + 1;
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

    % Trimmed RMSE: remove top 5% outliers for robust estimation
    trimPct = 95;
    velRMSE(m) = sqrt(mean(vErr2(vErr2 <= prctile(vErr2, trimPct))));
    rangeRMSE(m) = sqrt(mean(rErr2(rErr2 <= prctile(rErr2, trimPct))));
    velRMSE_ofdm(m) = sqrt(mean(vErr2_o(vErr2_o <= prctile(vErr2_o, trimPct))));
    rangeRMSE_ofdm(m) = sqrt(mean(rErr2_o(rErr2_o <= prctile(rErr2_o, trimPct))));
    fprintf('  EbNo=%+3d: OTFS vel=%.2f/rng=%.3f  OFDM vel=%.2f/rng=%.3f\n', ...
        EbNo_nav(m), velRMSE(m), rangeRMSE(m), velRMSE_ofdm(m), rangeRMSE_ofdm(m));
end

% CRB
T_frame = ofdmSym * Ts_sym;
snr_lin = 10.^(snr_nav/10);
crb_vel = (1 ./ (2*pi*T_frame*sqrt(2*snr_lin))) * c_light / fc;
crb_range = (1 ./ (2*pi*Bw*sqrt(2*snr_lin))) * c_light;

figure('Name', 'Fig5: Nav Accuracy vs SNR', 'Position', [100 100 900 400]);

subplot(1,2,1);
semilogy(EbNo_nav, velRMSE_ofdm, 'g-o', 'LineWidth', 2, 'MarkerSize', 6); hold on;
semilogy(EbNo_nav, velRMSE, 'b-d', 'LineWidth', 2, 'MarkerSize', 7);
semilogy(EbNo_nav, crb_vel, 'r--', 'LineWidth', 1.5);
xlabel('E_b/N_0 (dB)'); ylabel('Velocity RMSE (m/s)');
title('(a) Velocity Estimation');
legend('OFDM Pilot (SFFT+Quinn)', 'OTFS DD Pilot (Fractional)', 'CRB', ...
    'Location', 'southwest');
grid on; set(gca, 'FontSize', 11);

subplot(1,2,2);
semilogy(EbNo_nav, rangeRMSE_ofdm, 'g-o', 'LineWidth', 2, 'MarkerSize', 6); hold on;
semilogy(EbNo_nav, rangeRMSE, 'b-d', 'LineWidth', 2, 'MarkerSize', 7);
semilogy(EbNo_nav, crb_range, 'r--', 'LineWidth', 1.5);
xlabel('E_b/N_0 (dB)'); ylabel('Range RMSE (m)');
title('(b) Range Estimation');
legend('OFDM Pilot (SFFT+Quinn)', 'OTFS DD Pilot (Fractional)', 'CRB', ...
    'Location', 'southwest');
grid on; set(gca, 'FontSize', 11);

saveas(gcf, 'fig5_nav_accuracy_snr.fig');
print(gcf, 'fig5_nav_accuracy_snr', '-dpng', '-r150');
fprintf('  Fig 5 done (%.1f s)\n\n', toc(tFig));

%% ========================================================================
%  Fig 6: Navigation Accuracy vs Elevation — OTFS-Pilot vs OFDM-Pilot
%  ========================================================================
fprintf('=== Fig 6: Navigation vs Elevation (OTFS vs OFDM) ===\n');
tFig = tic;

fixedEbNo6 = 15;
snr6 = fixedEbNo6 + 10*log10(codeRate*k) + 10*log10(numDC/numSC) ;
numTrials6 = 200;

velRMSE_elev_otfs  = zeros(numElev, 1);
rangeRMSE_elev_otfs = zeros(numElev, 1);
velRMSE_elev_ofdm  = zeros(numElev, 1);
rangeRMSE_elev_ofdm = zeros(numElev, 1);

for ei = 1:numElev
    ntnCfg_e = ntnConfig;
    ntnCfg_e.elevAngle = elevAngles(ei);

    vErr2_t = zeros(numTrials6, 1);  rErr2_t = zeros(numTrials6, 1);
    vErr2_o = zeros(numTrials6, 1);  rErr2_o = zeros(numTrials6, 1);

    parfor trial = 1:numTrials6
        [chTF6, chN] = multipathChannel(cpSize, scs, tfGridSize, velocity, ntnCfg_e);

        % ---- Add random residual propagation delay ----
        bulkDelayBins = 2 + 6*rand();
        bulkDelay_s   = bulkDelayBins * delta_tau_s;
        chN.pathDelays_s = chN.pathDelays_s + bulkDelay_s;

        sc_idx_6 = (0:numSC-1)';
        phaseRamp6 = exp(-2j*pi * sc_idx_6 * scs * bulkDelay_s);
        chTF6 = chTF6 .* phaseRamp6;

        % True LoS (always first tap by construction in multipathChannel)
        trueVel   = chN.pathDopplers_Hz(1) * c_light / fc;
        trueRange = chN.pathDelays_s(1) * c_light;

        % ---- OTFS DD-Pilot sensing ----
        [dG, ~, ~, ~, pI] = pilotPatternDD(qamTx_pilot, ddGridSize(1), ofdmSym, pilotCfg);
        pI.maxDelayBins = 10;   % Cover bulk delay (2-8 bins) + multipath
        pI.maxDopplerBins = pilotCfg.maxDopplerBins;
        rF = applyChannelDD(dG, chN, scs, cpSize, fd_hz, true);
        sp = mean(abs(rF(:)).^2);
        nV = sp * 10^(-snr6/10);
        rF = rF + sqrt(nV/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));

        [~,~,~,nE] = ddChannelEstimate(rF, pI);

        estVel_t   = (nE.losFracDoppler * delta_nu_hz + fd_hz) * c_light / fc;
        estRange_t = nE.losFracDelay * delta_tau_s * c_light;

        vErr2_t(trial) = (estVel_t - trueVel)^2;
        rErr2_t(trial) = (estRange_t - trueRange)^2;

        % ---- OFDM TF-Pilot sensing (SFFT) ----
        sp_o = mean(abs(chTF6(:)).^2);
        nV_o = sp_o * 10^(-snr6/10);
        chTF_noisy = chTF6 + sqrt(nV_o/2) * ...
            (randn(size(chTF6)) + 1j*randn(size(chTF6)));

        % Doppler pre-compensation in TF domain: remove bulk fd
        m_precomp6 = exp(-1j*2*pi*fd_hz*(1:ofdmSym)*Ts_sym);
        chTF_precomp6 = chTF_noisy .* m_precomp6;

        % Symplectic FFT: TF → delay-Doppler response
        ddResp = fftshift(fft(ifft(chTF_precomp6, [], 1), [], 2));

        % Find LoS peak
        [~, peakIdx] = max(abs(ddResp(:)));
        [peakRow, peakCol] = ind2sub(size(ddResp), peakIdx);

        % Delay & Doppler estimation with Quinn fractional refinement
        delayCentre = floor(numSC/2) + 1;
        estDelay_o_s = (peakRow - delayCentre) * delta_tau_s;

        dopplerCentre = floor(ofdmSym/2) + 1;
        estResidualDoppler_o = (peakCol - dopplerCentre) * delta_nu_hz;

        % Doppler Quinn
        if peakCol > 1 && peakCol < ofdmSym
            Xm1 = ddResp(peakRow, peakCol-1);
            X0  = ddResp(peakRow, peakCol);
            Xp1 = ddResp(peakRow, peakCol+1);
            ap = real(Xp1/X0); am = real(Xm1/X0);
            dp = -ap/(1-ap); dm = am/(1-am);
            if abs(dp) < abs(dm), dk = dp; else, dk = dm; end
            dk = max(-0.5, min(0.5, dk));
            estResidualDoppler_o = (peakCol - dopplerCentre + dk) * delta_nu_hz;
        end
        % Delay Quinn
        if peakRow > 1 && peakRow < numSC
            Xm1 = ddResp(peakRow-1, peakCol);
            X0  = ddResp(peakRow, peakCol);
            Xp1 = ddResp(peakRow+1, peakCol);
            ap = real(Xp1/X0); am = real(Xm1/X0);
            dp = -ap/(1-ap); dm = am/(1-am);
            if abs(dp) < abs(dm), dl = dp; else, dl = dm; end
            dl = max(-0.5, min(0.5, dl));
            estDelay_o_s = (peakRow - delayCentre + dl) * delta_tau_s;
        end

        % Add back pre-compensated bulk Doppler
        estDoppler_o_Hz = estResidualDoppler_o + fd_hz;
        estVel_o = estDoppler_o_Hz * c_light / fc;
        estRange_o = abs(estDelay_o_s) * c_light;

        vErr2_o(trial) = (estVel_o - trueVel)^2;
        rErr2_o(trial) = (estRange_o - trueRange)^2;
    end

    % Trimmed RMSE: remove top 5% outliers for robust estimation
    trimPct = 95;
    velRMSE_elev_otfs(ei)  = sqrt(mean(vErr2_t(vErr2_t <= prctile(vErr2_t, trimPct))));
    rangeRMSE_elev_otfs(ei) = sqrt(mean(rErr2_t(rErr2_t <= prctile(rErr2_t, trimPct))));
    velRMSE_elev_ofdm(ei)  = sqrt(mean(vErr2_o(vErr2_o <= prctile(vErr2_o, trimPct))));
    rangeRMSE_elev_ofdm(ei) = sqrt(mean(rErr2_o(rErr2_o <= prctile(rErr2_o, trimPct))));

    fprintf('  Elev=%2d°: OTFS vel=%.2f/rng=%.3f  OFDM vel=%.1f/rng=%.3f\n', ...
        elevAngles(ei), velRMSE_elev_otfs(ei), rangeRMSE_elev_otfs(ei), ...
        velRMSE_elev_ofdm(ei), rangeRMSE_elev_ofdm(ei));
end

% ---- Plot Fig 6: 2x2 layout ----
figure('Name', 'Fig6: Nav vs Elevation', 'Position', [100 100 900 700]);

% (a) Velocity RMSE vs Elevation
subplot(2,2,1);
semilogy(elevAngles, velRMSE_elev_ofdm, 'g-o', 'LineWidth', 2, 'MarkerSize', 7); hold on;
semilogy(elevAngles, velRMSE_elev_otfs, 'b-d', 'LineWidth', 2, 'MarkerSize', 7);
xlabel('Elevation Angle (\circ)'); ylabel('Velocity RMSE (m/s)');
title(sprintf('(a) Velocity RMSE (E_b/N_0=%d dB)', fixedEbNo6));
legend('OFDM (SFFT+Quinn)', 'OTFS DD Pilot', 'Location', 'northeast');
grid on; xticks(elevAngles); set(gca, 'FontSize', 11);

% (b) Range RMSE vs Elevation
subplot(2,2,2);
semilogy(elevAngles, rangeRMSE_elev_ofdm, 'g-o', 'LineWidth', 2, 'MarkerSize', 7); hold on;
semilogy(elevAngles, rangeRMSE_elev_otfs, 'b-d', 'LineWidth', 2, 'MarkerSize', 7);
xlabel('Elevation Angle (\circ)'); ylabel('Range RMSE (m)');
title(sprintf('(b) Range RMSE (E_b/N_0=%d dB)', fixedEbNo6));
legend('OFDM (SFFT+Quinn)', 'OTFS DD Pilot', 'Location', 'northeast');
grid on; xticks(elevAngles); set(gca, 'FontSize', 11);

% (c) Velocity improvement ratio
subplot(2,2,3);
velGain = velRMSE_elev_ofdm ./ velRMSE_elev_otfs;
bar(elevAngles, velGain, 0.5, 'FaceColor', [0.20 0.40 0.85], 'EdgeColor', [0.12 0.24 0.55]);
hold on;
yline(1, 'k--', 'LineWidth', 1);
xlabel('Elevation Angle (\circ)'); ylabel('OFDM / OTFS Ratio');
title('(c) Velocity: OTFS Advantage Factor');
grid on; xticks(elevAngles); set(gca, 'FontSize', 11);
for j = 1:numElev
    text(elevAngles(j), velGain(j)+max(velGain)*0.03, ...
        sprintf('%.0f\\times', velGain(j)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, ...
        'FontWeight', 'bold', 'Color', [0.12 0.24 0.55]);
end

% (d) Range improvement ratio
subplot(2,2,4);
rngGain = rangeRMSE_elev_ofdm ./ rangeRMSE_elev_otfs;
bar(elevAngles, rngGain, 0.5, 'FaceColor', [0.85 0.33 0.10], 'EdgeColor', [0.55 0.20 0.06]);
hold on;
yline(1, 'k--', 'LineWidth', 1);
xlabel('Elevation Angle (\circ)'); ylabel('OFDM / OTFS Ratio');
title('(d) Range: OTFS Advantage Factor');
grid on; xticks(elevAngles); set(gca, 'FontSize', 11);
for j = 1:numElev
    text(elevAngles(j), rngGain(j)+max(rngGain)*0.03, ...
        sprintf('%.1f\\times', rngGain(j)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, ...
        'FontWeight', 'bold', 'Color', [0.55 0.20 0.06]);
end

sgtitle(sprintf('Navigation: OTFS-Pilot vs OFDM-Pilot (E_b/N_0=%d dB, NTN-TDL-C)', fixedEbNo6), ...
    'FontSize', 13, 'FontWeight', 'bold');

saveas(gcf, 'fig6_nav_vs_elevation.fig');
print(gcf, 'fig6_nav_vs_elevation', '-dpng', '-r150');
fprintf('  Fig 6 done (%.1f s)\n\n', toc(tFig));

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
fprintf('\n7 figures generated. Total time: %.1f s (%.1f min)\n', toc(tTotal), toc(tTotal)/60);
