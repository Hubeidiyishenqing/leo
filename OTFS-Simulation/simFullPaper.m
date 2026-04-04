%--------------------------------------------------------------------------
%
%  Comprehensive Simulation for Paper:
%    "DD-Domain Pilot Design for OTFS-Based ISAC in LEO NTN"
%
%  Generates Figures:
%    Fig 1: Four-scheme BER vs Eb/No (OFDM, OFDM-SIC, OTFS-Ideal, OTFS-Pilot)
%    Fig 4: Pilot Overhead vs Pre-compensation Ratio (analytical)
%    Fig 5: Velocity/Range RMSE vs SNR — with/without Quinn (ablation)
%    Fig 6: Sensing RMSE vs Ephemeris Error (robustness)
%    Fig 7: DD Grid Pilot Pattern Visualization
%    Fig 8: Sensing RMSE vs Elevation Angle / K-factor
%
%  Optimized: Split Fig1 into 2 parfor passes for memory efficiency
%--------------------------------------------------------------------------

clear; close all;
rng(42);
tTotal = tic;

% Initialize parallel pool — use most logical cores (split parfor keeps memory low)
numWorkers = max(2, feature('numcores') - 2);  % 12逻辑核 → 10 workers
if isempty(gcp('nocreate'))
    try
        parpool('Processes', numWorkers);
        fprintf('Parallel pool started with %d workers.\n', gcp().NumWorkers);
    catch ME
        warning('Parallel pool failed: %s\nRunning sequentially.', ME.message);
    end
else
    fprintf('Using existing parallel pool (%d workers).\n', gcp().NumWorkers);
end

%% ========================================================================
%  System Parameters
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

% QAM constellation
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
subframeBits = numDC * ofdmSym * k;
numCW = floor(subframeBits / noCodedbits);
if numCW < 1, numCW = 1; end
codedBitsTotal = numCW * noCodedbits;
padBits = max(0, subframeBits - codedBitsTotal);

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
%  OTFS-Pilot Specific Data (DD-domain)
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
%  Fig 1: Four-Scheme BER vs Eb/No  (Split into 2 passes for memory)
%  ========================================================================
fprintf('\n=== Fig 1: Four-Scheme BER vs Eb/No ===\n');
tFig = tic;

EbNo_fig1 = (0:2:24)';
numEbNo1 = length(EbNo_fig1);
minTrials = 50;          % Minimum MC trials per SNR point
maxTrials = 200;         % Maximum MC trials (for low BER)
minErrors = 100;         % Target minimum error events for reliability

% SNR formulas
snr_ofdm = EbNo_fig1 + 10*log10(codeRate*k) + 10*log10(numDC/numSC);
snr_otfs = snr_ofdm;
snr_pilot = EbNo_fig1 + 10*log10(codeRate*k) + 10*log10(numDataPosPilot/(numSC*ofdmSym));

% Uncoded BER
ber1_ofdm = zeros(numEbNo1, 1);
ber1_ofdm_sic = zeros(numEbNo1, 1);
ber1_otfs = zeros(numEbNo1, 1);
ber1_pilot = zeros(numEbNo1, 1);
% Coded BER (LDPC)
ber1_c_ofdm = zeros(numEbNo1, 1);
ber1_c_ofdm_sic = zeros(numEbNo1, 1);
ber1_c_otfs = zeros(numEbNo1, 1);
ber1_c_pilot = zeros(numEbNo1, 1);
trialCount = zeros(numEbNo1, 1);

for m = 1:numEbNo1
    snr_o_m = snr_ofdm(m);  snr_t_m = snr_otfs(m);  snr_p_m = snr_pilot(m);
    nBitsTx = length(txBits);  nBitsInfo = length(infoBits_in);
    nBitsPilot = length(txBits_pilot);  nBitsInfoP = length(infoBits_pilot);

    % ===== PASS A: OFDM + OFDM-SIC (adaptive MC) =====
    ber_o = 0; ber_os = 0; ber_co = 0; ber_cos = 0;
    totalTrials_A = 0;

    while totalTrials_A < maxTrials
        batchSize = min(minTrials, maxTrials - totalTrials_A);
        b_o=0; b_os=0; b_co=0; b_cos=0;

        parfor trial = 1:batchSize
            locDec = comm.LDPCDecoder(parityCheck_matrix);
            locDec.MaximumIterationCount = maxLDPCIter;
            locDec2 = comm.LDPCDecoder(parityCheck_matrix);
            locDec2.MaximumIterationCount = maxLDPCIter;

            [~, chI] = multipathChannel(cpSize, scs, tfGridSize, velocity, ntnConfig);
            chTF_ofdm = buildTF_OFDMeq(chI, scs, cpSize, numSC, ofdmSym, fd_hz);

            % --- Uncoded OFDM ---
            fTF = applyChannelTF_ICI(guardbandTx, chI, scs, cpSize, fd_hz);
            sp = mean(abs(fTF(:)).^2);
            nV = sp / 10^(snr_o_m/10);
            ns = sqrt(nV/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
            rxTF_ofdm = fTF + ns;
            eq = conj(chTF_ofdm) ./ (abs(chTF_ofdm).^2 + nV) .* rxTF_ofdm;
            pR = eq; pR(numDC/2+1:numDC/2+11,:)=[]; pR(1,:)=[];
            rb = qamdemod(pR(:), ModOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
            b_o = b_o + sum(txBits ~= rb(1:nBitsTx)) / nBitsTx;

            % --- Uncoded OFDM-SIC ---
            eq_sic = applyOFDM_ICICancel(rxTF_ofdm, chI, scs, cpSize, fd_hz, numSC, ofdmSym, nV, 1);
            pR_sic = eq_sic; pR_sic(numDC/2+1:numDC/2+11,:)=[]; pR_sic(1,:)=[];
            rb_sic = qamdemod(pR_sic(:), ModOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
            b_os = b_os + sum(txBits ~= rb_sic(1:nBitsTx)) / nBitsTx;

            % --- Coded OFDM ---
            fTF_c = applyChannelTF_ICI(guardbandTx_coded, chI, scs, cpSize, fd_hz);
            sp_c = mean(abs(fTF_c(:)).^2);
            nV_c = sp_c / 10^(snr_o_m/10);
            ns_c = sqrt(nV_c/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
            rxTF_ofdm_c = fTF_c + ns_c;

            eq_c = conj(chTF_ofdm) ./ (abs(chTF_ofdm).^2 + nV_c) .* rxTF_ofdm_c;
            pR_c = eq_c; pR_c(numDC/2+1:numDC/2+11,:)=[]; pR_c(1,:)=[];
            hardDec_o = qammod(qamdemod(pR_c(:), ModOrder, 'UnitAveragePower', true), ModOrder, 'UnitAveragePower', true);
            nV_c_emp = max(mean(abs(pR_c(:) - hardDec_o).^2), 1e-10);
            llr = qamdemod(pR_c(:), ModOrder, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', nV_c_emp);
            deintLLR = randdeintrlv(llr, 4831);
            deintLLR = deintLLR(1:numCW*noCodedbits);
            deintLLR = max(min(deintLLR, 50), -50);
            decBits = decodeLDPC_helper(deintLLR, locDec, numCW, noCodedbits);
            b_co = b_co + sum(infoBits_in ~= decBits(1:nBitsInfo)) / nBitsInfo;

            % --- Coded OFDM-SIC ---
            eq_sic_c = applyOFDM_ICICancel(rxTF_ofdm_c, chI, scs, cpSize, fd_hz, numSC, ofdmSym, nV_c, 1);
            pR_sic_c = eq_sic_c; pR_sic_c(numDC/2+1:numDC/2+11,:)=[]; pR_sic_c(1,:)=[];
            hardDec_os = qammod(qamdemod(pR_sic_c(:), ModOrder, 'UnitAveragePower', true), ModOrder, 'UnitAveragePower', true);
            nV_sic_emp = max(mean(abs(pR_sic_c(:) - hardDec_os).^2), 1e-10);
            llr_sic = qamdemod(pR_sic_c(:), ModOrder, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', nV_sic_emp);
            deintLLR_sic = randdeintrlv(llr_sic, 4831);
            deintLLR_sic = deintLLR_sic(1:numCW*noCodedbits);
            deintLLR_sic = max(min(deintLLR_sic, 50), -50);
            decBits_sic = decodeLDPC_helper(deintLLR_sic, locDec2, numCW, noCodedbits);
            b_cos = b_cos + sum(infoBits_in ~= decBits_sic(1:nBitsInfo)) / nBitsInfo;
        end

        ber_o = ber_o + b_o; ber_os = ber_os + b_os;
        ber_co = ber_co + b_co; ber_cos = ber_cos + b_cos;
        totalTrials_A = totalTrials_A + batchSize;

        if totalTrials_A >= minTrials
            avgBER_A = (ber_co + ber_cos) / (2 * totalTrials_A);
            errCount_A = (ber_co + ber_cos) * nBitsInfo;
            if avgBER_A > 1e-3 || errCount_A >= minErrors
                break;
            end
        end
    end

    % ===== PASS B: OTFS-Ideal + OTFS-Pilot (adaptive MC) =====
    ber_t = 0; ber_p = 0; ber_ct = 0; ber_cp = 0;
    totalTrials_B = 0;

    while totalTrials_B < maxTrials
        batchSize = min(minTrials, maxTrials - totalTrials_B);
        b_t=0; b_p=0; b_ct=0; b_cp=0;

        parfor trial = 1:batchSize
            locDec3 = comm.LDPCDecoder(parityCheck_matrix);
            locDec3.MaximumIterationCount = maxLDPCIter;
            locDecP = comm.LDPCDecoder(parityCheck_matrix);
            locDecP.MaximumIterationCount = maxLDPCIter;

            [~, chI] = multipathChannel(cpSize, scs, tfGridSize, velocity, ntnConfig);
            chTF_eq = buildTF_AFC(chI, scs, cpSize, numSC, ofdmSym, fd_hz);

            % --- Uncoded OTFS-Ideal ---
            rxDD_t = applyChannelDD(guardbandTx, chI, scs, cpSize, fd_hz);
            fTF2 = ISFFT(rxDD_t);
            sp2 = mean(abs(fTF2(:)).^2);
            nV2 = sp2 / 10^(snr_t_m/10);
            ns2 = sqrt(nV2/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
            eq2 = conj(chTF_eq) ./ (abs(chTF_eq).^2 + nV2) .* (fTF2 + ns2);
            rDD = SFFT(eq2);
            pR2 = rDD; pR2(numDC/2+1:numDC/2+11,:)=[]; pR2(1,:)=[];
            rb2 = qamdemod(pR2(:), ModOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
            b_t = b_t + sum(txBits ~= rb2(1:nBitsTx)) / nBitsTx;

            % --- Coded OTFS-Ideal ---
            rxDD_tc = applyChannelDD(guardbandTx_coded, chI, scs, cpSize, fd_hz);
            fTF2_c = ISFFT(rxDD_tc);
            sp2_c = mean(abs(fTF2_c(:)).^2);
            nV2_c = sp2_c / 10^(snr_t_m/10);
            ns2_c = sqrt(nV2_c/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
            eq2_c = conj(chTF_eq) ./ (abs(chTF_eq).^2 + nV2_c) .* (fTF2_c + ns2_c);
            rDD_c = SFFT(eq2_c);
            pR2_c = rDD_c; pR2_c(numDC/2+1:numDC/2+11,:)=[]; pR2_c(1,:)=[];
            hardDec2 = qammod(qamdemod(pR2_c(:), ModOrder, 'UnitAveragePower', true), ModOrder, 'UnitAveragePower', true);
            nV2_emp = max(mean(abs(pR2_c(:) - hardDec2).^2), 1e-10);
            llr2 = qamdemod(pR2_c(:), ModOrder, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', nV2_emp);
            deintLLR2 = randdeintrlv(llr2, 4831);
            deintLLR2 = deintLLR2(1:numCW*noCodedbits);
            deintLLR2 = max(min(deintLLR2, 50), -50);
            decBits2 = decodeLDPC_helper(deintLLR2, locDec3, numCW, noCodedbits);
            b_ct = b_ct + sum(infoBits_in ~= decBits2(1:nBitsInfo)) / nBitsInfo;

            % --- Uncoded OTFS-Pilot (frac-aware cancel) ---
            [dG, dI, ~, ~, pI] = pilotPatternDD(qamTx_pilot, ddGridSize(1), ofdmSym, pilotCfg);
            pI.maxDelayBins = pilotCfg.maxDelayBins;
            pI.maxDopplerBins = pilotCfg.maxDopplerBins;
            rF = applyChannelDD(dG, chI, scs, cpSize, fd_hz);
            sp3 = mean(abs(rF(:)).^2);
            nV3 = sp3 / 10^(snr_p_m/10);
            ns3 = sqrt(nV3/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
            rF_noisy = rF + ns3;

            [~,~,~,nI] = ddChannelEstimate(rF_noisy, pI, true);
            pA = dG(pilotCfg.lp, pilotCfg.kp);
            cG = nI.pathGains / pA;
            estCh = struct('pathDelays_s', nI.pathDelays * delta_tau_s, ...
                'pathDopplers_Hz', nI.pathDopplers * delta_nu_hz + fd_hz, ...
                'pathGains', cG, 'numPaths', nI.numPathsDetected);
            chTF_est = buildTF_AFC(estCh, scs, cpSize, numSC, ofdmSym, fd_hz);
            rClean = rF_noisy - nI.pilotResponse;
            rClean_TF = ISFFT(rClean);
            eq_p = conj(chTF_est) ./ (abs(chTF_est).^2 + nV3) .* rClean_TF;
            rDD_eq = SFFT(eq_p);
            eqS = rDD_eq(dI);
            rb_p = qamdemod(eqS, ModOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
            b_p = b_p + sum(txBits_pilot ~= rb_p(1:nBitsPilot)) / nBitsPilot;

            % --- Coded OTFS-Pilot ---
            [dG_c, dI_c, ~, ~, pI_c] = pilotPatternDD(qamTx_pilot_coded, ddGridSize(1), ofdmSym, pilotCfg);
            pI_c.maxDelayBins = pilotCfg.maxDelayBins;
            pI_c.maxDopplerBins = pilotCfg.maxDopplerBins;
            rF_c = applyChannelDD(dG_c, chI, scs, cpSize, fd_hz);
            sp3_c = mean(abs(rF_c(:)).^2);
            nV3_c = sp3_c / 10^(snr_p_m/10);
            ns3_c = sqrt(nV3_c/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));
            rF_c_noisy = rF_c + ns3_c;

            [~,~,~,nI_c] = ddChannelEstimate(rF_c_noisy, pI_c, true);
            pA_c = dG_c(pilotCfg.lp, pilotCfg.kp);
            cG_c = nI_c.pathGains / pA_c;
            estCh_c = struct('pathDelays_s', nI_c.pathDelays * delta_tau_s, ...
                'pathDopplers_Hz', nI_c.pathDopplers * delta_nu_hz + fd_hz, ...
                'pathGains', cG_c, 'numPaths', nI_c.numPathsDetected);
            chTF_est_c = buildTF_AFC(estCh_c, scs, cpSize, numSC, ofdmSym, fd_hz);
            rClean_c = rF_c_noisy - nI_c.pilotResponse;
            rClean_c_TF = ISFFT(rClean_c);
            eq_pc = conj(chTF_est_c) ./ (abs(chTF_est_c).^2 + nV3_c) .* rClean_c_TF;
            rDD_eq_c = SFFT(eq_pc);
            eqS_c = rDD_eq_c(dI_c);
            hardDec_p = qammod(qamdemod(eqS_c, ModOrder, 'UnitAveragePower', true), ModOrder, 'UnitAveragePower', true);
            nV3_c_emp = max(mean(abs(eqS_c - hardDec_p).^2), 1e-10);
            llr_p = qamdemod(eqS_c, ModOrder, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', nV3_c_emp);
            deintLLR_p = randdeintrlv(llr_p, 4831);
            deintLLR_p = deintLLR_p(1:numCW_pilot*noCodedbits);
            deintLLR_p = max(min(deintLLR_p, 50), -50);
            decBits_p = decodeLDPC_helper(deintLLR_p, locDecP, numCW_pilot, noCodedbits);
            b_cp = b_cp + sum(infoBits_pilot ~= decBits_p(1:nBitsInfoP)) / nBitsInfoP;
        end

        ber_t = ber_t + b_t; ber_p = ber_p + b_p;
        ber_ct = ber_ct + b_ct; ber_cp = ber_cp + b_cp;
        totalTrials_B = totalTrials_B + batchSize;

        if totalTrials_B >= minTrials
            avgBER_B = (ber_ct + ber_cp) / (2 * totalTrials_B);
            errCount_B = (ber_ct + ber_cp) * nBitsInfo;
            if avgBER_B > 1e-3 || errCount_B >= minErrors
                break;
            end
        end
    end

    % Use the larger trial count for normalization
    numTrials1 = max(totalTrials_A, totalTrials_B);
    ber1_ofdm(m) = ber_o / totalTrials_A;
    ber1_ofdm_sic(m) = ber_os / totalTrials_A;
    ber1_otfs(m) = ber_t / totalTrials_B;
    ber1_pilot(m) = ber_p / totalTrials_B;
    ber1_c_ofdm(m) = ber_co / totalTrials_A;
    ber1_c_ofdm_sic(m) = ber_cos / totalTrials_A;
    ber1_c_otfs(m) = ber_ct / totalTrials_B;
    ber1_c_pilot(m) = ber_cp / totalTrials_B;
    trialCount(m) = numTrials1;

    fprintf('  EbNo=%+3d dB [A=%3d,B=%3d]: OFDM=%.4f/%.4f  SIC=%.4f/%.4f  OTFS=%.4f/%.4f  Pilot=%.4f/%.4f\n', ...
        EbNo_fig1(m), totalTrials_A, totalTrials_B, ...
        ber1_ofdm(m), ber1_c_ofdm(m), ...
        ber1_ofdm_sic(m), ber1_c_ofdm_sic(m), ...
        ber1_otfs(m), ber1_c_otfs(m), ...
        ber1_pilot(m), ber1_c_pilot(m));
end

% Plot Fig 1: 8 curves (4 uncoded solid + 4 coded dashed)
figure('Name', 'Fig1: Four-Scheme BER', 'Position', [100 100 800 620]);
% Uncoded (solid)
semilogy(EbNo_fig1, max(ber1_ofdm, 1e-6), 'g-o', 'LineWidth', 2.5, 'MarkerSize', 9); hold on;
semilogy(EbNo_fig1, max(ber1_ofdm_sic, 1e-6), 'm-^', 'LineWidth', 2.5, 'MarkerSize', 9);
semilogy(EbNo_fig1, max(ber1_otfs, 1e-6), 'r-s', 'LineWidth', 2.5, 'MarkerSize', 9);
semilogy(EbNo_fig1, max(ber1_pilot, 1e-6), 'b-d', 'LineWidth', 2.5, 'MarkerSize', 9);
% Coded (dashed)
semilogy(EbNo_fig1, max(ber1_c_ofdm, 1e-6), 'g--o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
semilogy(EbNo_fig1, max(ber1_c_ofdm_sic, 1e-6), 'm--^', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'm');
semilogy(EbNo_fig1, max(ber1_c_otfs, 1e-6), 'r--s', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'r');
semilogy(EbNo_fig1, max(ber1_c_pilot, 1e-6), 'b--d', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
xlabel('$E_b/N_0$ (dB)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('BER', 'FontSize', 14);
legend('OFDM', 'OFDM-SIC', 'OTFS (ideal CSI)', 'OTFS-Pilot', ...
    'C-OFDM', 'C-OFDM-SIC', 'C-OTFS', 'C-OTFS-Pilot', ...
    'Location', 'southwest', 'FontSize', 10, 'NumColumns', 2);
grid on; ylim([1e-6, 1]);
set(gca, 'FontSize', 13, 'LineWidth', 1);

saveas(gcf, 'fig1_three_scheme_ber.fig');
print(gcf, 'fig1_three_scheme_ber', '-dpng', '-r300');
fprintf('  Fig 1 done (%.1f s)\n\n', toc(tFig));

% Print statistics table
fprintf('  --- BER Statistics ---\n');
fprintf('  %6s  %6s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n', ...
    'EbNo', 'Trials', 'OFDM-U', 'OFDM-SIC-U', 'OTFS-U', 'Pilot-U', ...
    'OFDM-C', 'OFDM-SIC-C', 'OTFS-C', 'Pilot-C');
for m = 1:numEbNo1
    fprintf('  %+4d dB  %4d    %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e\n', ...
        EbNo_fig1(m), trialCount(m), ...
        ber1_ofdm(m), ber1_ofdm_sic(m), ber1_otfs(m), ber1_pilot(m), ...
        ber1_c_ofdm(m), ber1_c_ofdm_sic(m), ber1_c_otfs(m), ber1_c_pilot(m));
end
fprintf('\n');


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
isFeasible = false(numRatios, 1);

for idx = 1:numRatios
    ratio = precompRatios(idx);
    res_hz = (1 - ratio) * fd_hz + 0.05 * fd_hz;
    res_bins = ceil(res_hz / delta_nu_hz);
    kG = 2 * res_bins + 2;
    kGuardVals(idx) = kG;
    isFeasible(idx) = (kG < floor(ofdmSym/2));
    overheadPct(idx) = 100 * (2*kG+1) * (2*lGuardFixed+1) / (numSC * ofdmSym);
end

feasMask = isFeasible;
dataCap = max(0, 1 - overheadPct / 100);

feasThresh = NaN;
for idx = 1:numRatios
    if feasMask(idx), feasThresh = precompRatios(idx)*100; break; end
end

figure('Name', 'Fig4: Pilot Overhead', 'Position', [100 100 960 420]);

subplot(1,2,1);
plot(precompRatios*100, overheadPct, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 6);
hold on;
if ~isnan(feasThresh)
    xline(feasThresh, 'r--', sprintf('Feasibility (%.0f%%)', feasThresh), ...
        'LineWidth', 1.5, 'FontSize', 12, ...
        'LabelOrientation', 'aligned', 'LabelVerticalAlignment', 'bottom');
    patch([0 feasThresh feasThresh 0], [0 0 100 100], 'r', 'FaceAlpha', 0.08, 'EdgeColor', 'none');
end
xlabel('Doppler Pre-comp Ratio (%)', 'FontSize', 13);
ylabel('Overhead (%)', 'FontSize', 13);
title('(a) Pilot + Guard Overhead', 'FontSize', 14);
grid on;
yMax = max(100, 1.05 * max(overheadPct));
ylim([0 yMax]);
set(gca, 'FontSize', 13, 'LineWidth', 1);

subplot(1,2,2);
yyaxis left;
barIdx = [];
for t = [0 20 50 70 80 90 100]
    [~,bi] = min(abs(precompRatios*100-t)); barIdx = [barIdx, bi];
end
bar(precompRatios(barIdx)*100, kGuardVals(barIdx), 0.5);
ylabel('k_{Guard} (bins)', 'FontSize', 13);
hold on;
yline(floor(ofdmSym/2), 'k--', 'M/2 boundary', 'LineWidth', 1.2, ...
    'FontSize', 11, 'LabelVerticalAlignment', 'bottom');
yyaxis right;
plot(precompRatios*100, dataCap*100, 'r-s', 'LineWidth', 2.5, 'MarkerSize', 6);
ylabel('Data Capacity (%)', 'FontSize', 13);
xlabel('Pre-comp Ratio (%)', 'FontSize', 13);
title('(b) Guard Size & Capacity', 'FontSize', 14);
legend('k_{Guard}', 'M/2 boundary', 'Capacity', 'Location', 'east', 'FontSize', 12);
grid on; set(gca, 'FontSize', 13, 'LineWidth', 1);

saveas(gcf, 'fig4_overhead_vs_precomp.fig');
print(gcf, 'fig4_overhead_vs_precomp', '-dpng', '-r300');
fprintf('  Fig 4 done (%.1f s)\n\n', toc(tFig));

%% ========================================================================
%  Fig 5: Navigation RMSE — With vs Without Quinn (Ablation Study)
%  ========================================================================
fprintf('=== Fig 5: Navigation Accuracy — Quinn Ablation ===\n');
tFig = tic;

EbNo_nav = (0:2:24)';
numNavEbNo = length(EbNo_nav);
numNavTrials = 300;   % Reduced: median-based RMSE converges fast
snr_nav = EbNo_nav + 10*log10(codeRate*k) + 10*log10(numDC/numSC);

% With Quinn (fractional interpolation)
velRMSE = zeros(numNavEbNo, 1);
rangeRMSE = zeros(numNavEbNo, 1);
% Without Quinn (integer-only — ablation baseline)
velRMSE_noQ = zeros(numNavEbNo, 1);
rangeRMSE_noQ = zeros(numNavEbNo, 1);

for m = 1:numNavEbNo
    vErr2 = zeros(numNavTrials, 1);
    rErr2 = zeros(numNavTrials, 1);
    vErr2_noQ = zeros(numNavTrials, 1);
    rErr2_noQ = zeros(numNavTrials, 1);

    snr_nav_m = snr_nav(m);
    parfor trial = 1:numNavTrials
        [~, chN] = multipathChannel(cpSize, scs, tfGridSize, velocity, ntnConfig);

        bulkDelayBins = 2 + 6*rand();
        bulkDelay_s   = bulkDelayBins * delta_tau_s;
        chN.pathDelays_s = chN.pathDelays_s + bulkDelay_s;

        trueVel = chN.pathDopplers_Hz(1) * c_light / fc;
        trueRange = chN.pathDelays_s(1) * c_light;

        % OTFS-Pilot sensing
        [dG, ~, ~, ~, pI] = pilotPatternDD(qamTx_pilot, ddGridSize(1), ofdmSym, pilotCfg);
        pI.maxDelayBins = 10;
        pI.maxDopplerBins = pilotCfg.maxDopplerBins;
        rF = applyChannelDD(dG, chN, scs, cpSize, fd_hz, true);

        sp = mean(abs(rF(:)).^2);
        nV = sp * 10^(-snr_nav_m/10);
        rF = rF + sqrt(nV/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));

        [~,~,~,nE] = ddChannelEstimate(rF, pI, false);

        % WITH Quinn: fractional estimates
        estDoppler_Hz = nE.losFracDoppler * delta_nu_hz + fd_hz;
        estDelay_s = nE.losFracDelay * delta_tau_s;
        estVel = estDoppler_Hz * c_light / fc;
        estRange = estDelay_s * c_light;

        vErr2(trial) = (estVel - trueVel)^2;
        rErr2(trial) = (estRange - trueRange)^2;

        % WITHOUT Quinn: integer-only estimates (ablation)
        losIdx_nq = 1;
        maxAmp_nq = max(abs(nE.pathGains));
        strongMask_nq = abs(nE.pathGains) > maxAmp_nq * 10^(-6/20);
        strongIdx_nq = find(strongMask_nq);
        if ~isempty(strongIdx_nq)
            [~, minPos_nq] = min(abs(nE.pathDelays(strongIdx_nq)));
            losIdx_nq = strongIdx_nq(minPos_nq);
        end
        estDoppler_noQ = nE.pathDopplers(losIdx_nq) * delta_nu_hz + fd_hz;
        estDelay_noQ = nE.pathDelays(losIdx_nq) * delta_tau_s;
        estVel_noQ = estDoppler_noQ * c_light / fc;
        estRange_noQ = estDelay_noQ * c_light;

        vErr2_noQ(trial) = (estVel_noQ - trueVel)^2;
        rErr2_noQ(trial) = (estRange_noQ - trueRange)^2;
    end

    velRMSE(m) = sqrt(median(vErr2));
    rangeRMSE(m) = sqrt(median(rErr2));
    velRMSE_noQ(m) = sqrt(median(vErr2_noQ));
    rangeRMSE_noQ(m) = sqrt(median(rErr2_noQ));
    fprintf('  EbNo=%+3d dB: vel=%.2f/%.2f m/s, rng=%.3f/%.3f m  (Quinn/NoQuinn)\n', ...
        EbNo_nav(m), velRMSE(m), velRMSE_noQ(m), rangeRMSE(m), rangeRMSE_noQ(m));
end

% CRB
snr_lin = 10.^(snr_nav/10);
pilotBoost_lin = 10^(pilotCfg.pilotBoostdB/10);
crb_vel = (c_light/fc) * sqrt(12) ./ ...
    (2*pi*Ts_sym * sqrt(pilotBoost_lin .* snr_lin * ofdmSym * (ofdmSym^2-1)));
crb_range = c_light * sqrt(12) ./ ...
    (2*pi*scs * sqrt(pilotBoost_lin .* snr_lin * numSC * (numSC^2-1)));

figure('Name', 'Fig5: Nav Accuracy — Quinn Ablation', 'Position', [100 100 960 420]);

subplot(1,2,1);
semilogy(EbNo_nav, velRMSE, 'b-d', 'LineWidth', 2.5, 'MarkerSize', 9); hold on;
semilogy(EbNo_nav, velRMSE_noQ, 'k-x', 'LineWidth', 2, 'MarkerSize', 9);
semilogy(EbNo_nav, crb_vel, 'r--', 'LineWidth', 2);
xlabel('$E_b/N_0$ (dB)', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Velocity RMSE (m/s)', 'FontSize', 13);
title('(a) Velocity Estimation', 'FontSize', 14);
legend('With Quinn', 'Integer-only', 'CRB', 'Location', 'southwest', 'FontSize', 12);
grid on; set(gca, 'FontSize', 13, 'LineWidth', 1);

subplot(1,2,2);
semilogy(EbNo_nav, rangeRMSE, 'b-d', 'LineWidth', 2.5, 'MarkerSize', 9); hold on;
semilogy(EbNo_nav, rangeRMSE_noQ, 'k-x', 'LineWidth', 2, 'MarkerSize', 9);
semilogy(EbNo_nav, crb_range, 'r--', 'LineWidth', 2);
xlabel('$E_b/N_0$ (dB)', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Range RMSE (m)', 'FontSize', 13);
title('(b) Range Estimation', 'FontSize', 14);
legend('With Quinn', 'Integer-only', 'CRB', 'Location', 'southwest', 'FontSize', 12);
grid on; set(gca, 'FontSize', 13, 'LineWidth', 1);

saveas(gcf, 'fig5_nav_accuracy_snr.fig');
print(gcf, 'fig5_nav_accuracy_snr', '-dpng', '-r300');
fprintf('  Fig 5 done (%.1f s)\n\n', toc(tFig));


%% ========================================================================
%  Fig 6: Sensing RMSE vs Ephemeris Error (Robustness Analysis)
%  ========================================================================
fprintf('=== Fig 6: Sensing vs Ephemeris Error ===\n');
tFig = tic;

% Ephemeris error as fraction of fd
ephErrors = [0.0001, 0.001, 0.005, 0.01, 0.02, 0.05];  % 0.01% to 5%
numEph = length(ephErrors);
fixedEbNo = 20;  % dB
snr_eph = fixedEbNo + 10*log10(codeRate*k) + 10*log10(numDC/numSC);
numEphTrials = 1000;  % Need enough trials for monotonic trend

velRMSE_eph = zeros(numEph, 1);
rangeRMSE_eph = zeros(numEph, 1);

for e = 1:numEph
    vErr2 = zeros(numEphTrials, 1);
    rErr2 = zeros(numEphTrials, 1);

    eph_err = ephErrors(e);
    fd_precomp_eph = fd_hz * (1 - eph_err);

    res_hz_eph = eph_err * fd_hz + 0.05 * fd_hz;
    res_bins_eph = ceil(res_hz_eph / delta_nu_hz);
    pilotCfg_eph = pilotCfg;
    pilotCfg_eph.kGuard = min(2 * res_bins_eph + 2, floor(ofdmSym/2) - 1);
    pilotCfg_eph.maxDopplerBins = res_bins_eph + 1;

    parfor trial = 1:numEphTrials
        [~, chN] = multipathChannel(cpSize, scs, tfGridSize, velocity, ntnConfig);
        bulkDelayBins = 2 + 6*rand();
        chN.pathDelays_s = chN.pathDelays_s + bulkDelayBins * delta_tau_s;

        trueVel = chN.pathDopplers_Hz(1) * c_light / fc;
        trueRange = chN.pathDelays_s(1) * c_light;

        [dG, ~, ~, ~, pI] = pilotPatternDD(qamTx_pilot, ddGridSize(1), ofdmSym, pilotCfg_eph);
        pI.maxDelayBins = 10;
        pI.maxDopplerBins = pilotCfg_eph.maxDopplerBins;
        rF = applyChannelDD(dG, chN, scs, cpSize, fd_precomp_eph, true);

        sp = mean(abs(rF(:)).^2);
        nV = sp * 10^(-snr_eph/10);
        rF = rF + sqrt(nV/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));

        [~,~,~,nE] = ddChannelEstimate(rF, pI, false);

        estDoppler_Hz = nE.losFracDoppler * delta_nu_hz + fd_precomp_eph;
        estDelay_s = nE.losFracDelay * delta_tau_s;
        estVel = estDoppler_Hz * c_light / fc;
        estRange = estDelay_s * c_light;

        vErr2(trial) = (estVel - trueVel)^2;
        rErr2(trial) = (estRange - trueRange)^2;
    end

    velRMSE_eph(e) = sqrt(median(vErr2));
    rangeRMSE_eph(e) = sqrt(median(rErr2));
    fprintf('  Eph error=%.3f%%: vel=%.2f m/s, rng=%.3f m\n', ...
        eph_err*100, velRMSE_eph(e), rangeRMSE_eph(e));
end

figure('Name', 'Fig6: Sensing vs Ephemeris Error', 'Position', [100 100 960 420]);

subplot(1,2,1);
semilogx(ephErrors*100, velRMSE_eph, 'b-d', 'LineWidth', 2.5, 'MarkerSize', 9);
xlabel('Ephemeris Error (%)', 'FontSize', 13);
ylabel('Velocity RMSE (m/s)', 'FontSize', 13);
title('(a) Velocity vs Ephemeris Error', 'FontSize', 14);
grid on; set(gca, 'FontSize', 13, 'LineWidth', 1);
hold on; xline(5, 'r--', '5% boundary', 'LineWidth', 1.5, 'FontSize', 11, ...
    'LabelOrientation', 'aligned', 'LabelVerticalAlignment', 'bottom');

subplot(1,2,2);
semilogx(ephErrors*100, rangeRMSE_eph, 'b-d', 'LineWidth', 2.5, 'MarkerSize', 9);
xlabel('Ephemeris Error (%)', 'FontSize', 13);
ylabel('Range RMSE (m)', 'FontSize', 13);
title('(b) Range vs Ephemeris Error', 'FontSize', 14);
grid on; set(gca, 'FontSize', 13, 'LineWidth', 1);
hold on; xline(5, 'r--', '5% boundary', 'LineWidth', 1.5, 'FontSize', 11, ...
    'LabelOrientation', 'aligned', 'LabelVerticalAlignment', 'bottom');

saveas(gcf, 'fig6_eph_error_sweep.fig');
print(gcf, 'fig6_eph_error_sweep', '-dpng', '-r300');
fprintf('  Fig 6 done (%.1f s)\n\n', toc(tFig));


%% ========================================================================
%  Fig 7: DD Grid Pilot Pattern Visualization
%  ========================================================================
fprintf('=== Fig 7: DD Grid Visualization ===\n');
tFig = tic;

figure('Name', 'Fig7: DD Grid Comparison', 'Position', [100 100 1060 480]);

% (a) WITHOUT pre-comp
fullDop_bins = ceil(fd_hz / delta_nu_hz);
kG_noPre = 2 * fullDop_bins + 2;

subplot(1,2,1);
if kG_noPre >= floor(ofdmSym/2)
    imagesc(zeros(ddGridSize)); colormap(hot);
    text(ofdmSym/2, ddGridSize(1)/2, ...
        {'\bf INFEASIBLE', ...
         sprintf('k_{Guard}=%d > M/2=%d', kG_noPre, floor(ofdmSym/2)), ...
         'Guard exceeds grid boundary'}, ...
        'HorizontalAlignment', 'center', 'FontSize', 16, ...
        'Color', 'w', 'BackgroundColor', [0.8 0 0]);
    title(sprintf('(a) No Pre-comp (k_{Guard}=%d)', kG_noPre), 'FontSize', 14);
else
    pcNP = pilotCfg;
    pcNP.kGuard = kG_noPre;
    [~,dI_np,pI_np,gI_np,info_np] = pilotPatternDD(qamTx_pilot, ddGridSize(1), ofdmSym, pcNP);
    vis = zeros(ddGridSize); vis(dI_np)=0.3; vis(gI_np)=0; vis(pI_np)=1;
    imagesc(vis); colormap(hot);
    title(sprintf('(a) No Pre-comp: Overhead=%.1f%%', info_np.overheadPercent), 'FontSize', 14);
    hold on;
    rectangle('Position', [pcNP.kp-pcNP.kGuard-0.5, pcNP.lp-pcNP.lGuard-0.5, ...
        2*pcNP.kGuard+1, 2*pcNP.lGuard+1], 'EdgeColor','c','LineWidth',2.5,'LineStyle','--');
    hold off;
end
xlabel('Doppler Index', 'FontSize', 13); ylabel('Delay Index', 'FontSize', 13);
set(gca, 'FontSize', 13, 'LineWidth', 1);

% (b) WITH pre-comp
[~,dI_wp,pI_wp,gI_wp,info_wp] = pilotPatternDD(qamTx_pilot, ddGridSize(1), ofdmSym, pilotCfg);
subplot(1,2,2);
vis = zeros(ddGridSize); vis(dI_wp)=0.3; vis(gI_wp)=0; vis(pI_wp)=1;
imagesc(vis); colormap(hot);
title(sprintf('(b) With Pre-comp: Overhead=%.1f%%', info_wp.overheadPercent), 'FontSize', 14);
xlabel('Doppler Index', 'FontSize', 13); ylabel('Delay Index', 'FontSize', 13);
hold on;
rectangle('Position', [pilotCfg.kp-pilotCfg.kGuard-0.5, pilotCfg.lp-pilotCfg.lGuard-0.5, ...
    2*pilotCfg.kGuard+1, 2*pilotCfg.lGuard+1], 'EdgeColor','c','LineWidth',2.5,'LineStyle','--');
hold off;
set(gca, 'FontSize', 13, 'LineWidth', 1);

saveas(gcf, 'fig7_dd_grid_comparison.fig');
print(gcf, 'fig7_dd_grid_comparison', '-dpng', '-r300');
fprintf('  Fig 7 done (%.1f s)\n\n', toc(tFig));


%% ========================================================================
%  Fig 8: Sensing RMSE vs Elevation — Multi-Scenario (K-factor sweep)
%  DenseUrban K=4.4-6.8, Urban K=7-10, Suburban K=10-13 dB
%  Wide K-factor range shows clear RMSE trend
%  ========================================================================
fprintf('=== Fig 8: Sensing vs Elevation (Multi-Scenario) ===\n');
tFig = tic;

elevAngles = [20, 30, 40, 50, 60, 70, 80];
numElev = length(elevAngles);
scenarios = {'DenseUrban', 'Urban', 'Suburban'};
numScen = length(scenarios);
fixedEbNo_elev = 20;
snr_elev = fixedEbNo_elev + 10*log10(codeRate*k) + 10*log10(numDC/numSC);
numElevTrials = 1000;

velRMSE_elev = zeros(numElev, numScen);
rangeRMSE_elev = zeros(numElev, numScen);
kFactor_elev = zeros(numElev, numScen);

for s = 1:numScen
    fprintf('  --- Scenario: %s ---\n', scenarios{s});
    for e = 1:numElev
        ntnCfg_e = ntnConfig;
        ntnCfg_e.elevAngle = elevAngles(e);
        ntnCfg_e.scenario = scenarios{s};

        vErr2 = zeros(numElevTrials, 1);
        rErr2 = zeros(numElevTrials, 1);

        parfor trial = 1:numElevTrials
            [~, chN] = multipathChannel(cpSize, scs, tfGridSize, velocity, ntnCfg_e);
            bulkDelayBins = 2 + 6*rand();
            chN.pathDelays_s = chN.pathDelays_s + bulkDelayBins * delta_tau_s;

            trueVel = chN.pathDopplers_Hz(1) * c_light / fc;
            trueRange = chN.pathDelays_s(1) * c_light;

            [dG, ~, ~, ~, pI] = pilotPatternDD(qamTx_pilot, ddGridSize(1), ofdmSym, pilotCfg);
            pI.maxDelayBins = 10;
            pI.maxDopplerBins = pilotCfg.maxDopplerBins;
            rF = applyChannelDD(dG, chN, scs, cpSize, fd_hz, true);

            sp = mean(abs(rF(:)).^2);
            nV = sp * 10^(-snr_elev/10);
            rF = rF + sqrt(nV/2) * (randn(ddGridSize) + 1j*randn(ddGridSize));

            [~,~,~,nE] = ddChannelEstimate(rF, pI, false);

            estDoppler_Hz = nE.losFracDoppler * delta_nu_hz + fd_hz;
            estDelay_s = nE.losFracDelay * delta_tau_s;
            estVel = estDoppler_Hz * c_light / fc;
            estRange = estDelay_s * c_light;

            vErr2(trial) = (estVel - trueVel)^2;
            rErr2(trial) = (estRange - trueRange)^2;
        end

        % Trimmed mean (5%-95%) for outlier robustness
        vSorted = sort(vErr2); rSorted = sort(rErr2);
        trimN = round(0.05 * numElevTrials);
        trimIdx = (trimN+1):(numElevTrials-trimN);
        velRMSE_elev(e, s) = sqrt(mean(vSorted(trimIdx)));
        rangeRMSE_elev(e, s) = sqrt(mean(rSorted(trimIdx)));
        kFactor_elev(e, s) = getNtnKFactor_ext(elevAngles(e), scenarios{s});
        fprintf('    Elev=%d° (K=%.1f dB): vel=%.2f m/s, rng=%.3f m\n', ...
            elevAngles(e), kFactor_elev(e,s), velRMSE_elev(e,s), rangeRMSE_elev(e,s));
    end
end

% Colors/markers for each scenario
scenColors = {'r', [0 0.6 0], 'b'};   % DenseUrban=red, Urban=green, Suburban=blue
scenMarkers = {'s', '^', 'd'};

figure('Name', 'Fig8: Sensing vs Elevation (Multi-Scenario)', 'Position', [100 100 960 420]);

subplot(1,2,1);
hold on;
hLines = gobjects(numScen, 1);
for s = 1:numScen
    hLines(s) = plot(elevAngles, velRMSE_elev(:,s), '-', ...
        'Color', scenColors{s}, 'Marker', scenMarkers{s}, ...
        'LineWidth', 2.5, 'MarkerSize', 9);
end
xlabel('Elevation Angle ($^\circ$)', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Velocity RMSE (m/s)', 'FontSize', 13);
title('(a) Velocity Estimation', 'FontSize', 14);
legend(hLines, {'Dense Urban', 'Urban', 'Suburban'}, ...
    'Location', 'northeast', 'FontSize', 11);
grid on; set(gca, 'FontSize', 13, 'LineWidth', 1);

subplot(1,2,2);
hold on;
hLines2 = gobjects(numScen, 1);
for s = 1:numScen
    hLines2(s) = plot(elevAngles, rangeRMSE_elev(:,s), '-', ...
        'Color', scenColors{s}, 'Marker', scenMarkers{s}, ...
        'LineWidth', 2.5, 'MarkerSize', 9);
end
xlabel('Elevation Angle ($^\circ$)', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Range RMSE (m)', 'FontSize', 13);
title('(b) Range Estimation', 'FontSize', 14);
legend(hLines2, {'Dense Urban', 'Urban', 'Suburban'}, ...
    'Location', 'northeast', 'FontSize', 11);
grid on; set(gca, 'FontSize', 13, 'LineWidth', 1);

saveas(gcf, 'fig8_elev_sweep.fig');
print(gcf, 'fig8_elev_sweep', '-dpng', '-r300');
fprintf('  Fig 8 done (%.1f s)\n\n', toc(tFig));


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
    fprintf('    OFDM:       %.2e (coded: %.2e)\n', ber1_ofdm(idx20), ber1_c_ofdm(idx20));
    fprintf('    OFDM-SIC:   %.2e (coded: %.2e)\n', ber1_ofdm_sic(idx20), ber1_c_ofdm_sic(idx20));
    fprintf('    OTFS:       %.2e (coded: %.2e)\n', ber1_otfs(idx20), ber1_c_otfs(idx20));
    fprintf('    OTFS-Pilot: %.2e (coded: %.2e)\n', ber1_pilot(idx20), ber1_c_pilot(idx20));
end
fprintf('\n  --- Fig 5: Nav at Eb/No=20 dB ---\n');
idx20n = find(EbNo_nav==20, 1);
if ~isempty(idx20n)
    fprintf('    With Quinn:    vel=%.2f m/s, rng=%.3f m\n', velRMSE(idx20n), rangeRMSE(idx20n));
    fprintf('    Integer-only:  vel=%.2f m/s, rng=%.3f m\n', velRMSE_noQ(idx20n), rangeRMSE_noQ(idx20n));
end
fprintf('============================================================\n');
fprintf('\nAll figures generated. Total time: %.1f s (%.1f min)\n', toc(tTotal), toc(tTotal)/60);


%% ========================================================================
%  Local helper: K-factor lookup (duplicated from multipathChannel for
%  standalone access in the elevation sweep)
%  ========================================================================
function K_dB = getNtnKFactor_ext(elevAngle, scenario)
    elev_table = [10, 20, 30, 40, 50, 60, 70, 80, 90];
    K_Suburban_dB = [9.0, 10.0, 10.8, 11.4, 12.0, 12.4, 12.7, 12.9, 13.0];
    K_Urban_dB    = [7.0,  7.8,  8.3,  8.8,  9.2,  9.5,  9.7,  9.9, 10.0];
    K_Rural_dB    = [12.0, 13.2, 14.0, 14.6, 15.0, 15.4, 15.6, 15.8, 16.0];
    K_DenseUrban_dB = [4.4, 5.0, 5.4, 5.8, 6.2, 6.4, 6.6, 6.7, 6.8];
    switch scenario
        case 'Suburban', K_table = K_Suburban_dB;
        case 'Urban',    K_table = K_Urban_dB;
        case 'Rural',    K_table = K_Rural_dB;
        case 'DenseUrban', K_table = K_DenseUrban_dB;
        otherwise, K_table = K_Urban_dB;
    end
    elevAngle = max(10, min(90, elevAngle));
    K_dB = interp1(elev_table, K_table, elevAngle, 'pchip');
end


%% ========================================================================
%  Local helper: robust parallel pool startup
%  ========================================================================
function [poolObj, statusMsg] = ensureParallelPool(parallelCfg)
    poolObj = gcp('nocreate');
    if ~isempty(poolObj)
        statusMsg = sprintf('Reusing parallel pool: profile=%s, workers=%d.', ...
            poolObj.Cluster.Profile, poolObj.NumWorkers);
        return;
    end

    maxWorkers = max(1, floor(parallelCfg.maxWorkers));
    preferredProfile = parallelCfg.preferredProfile;

    try
        poolObj = parpool(preferredProfile, maxWorkers);
        statusMsg = sprintf('Parallel pool started: profile=%s, workers=%d.', ...
            preferredProfile, poolObj.NumWorkers);
        return;
    catch ME1
        statusMsg = sprintf(['Parallel pool startup failed for profile=%s ' ...
            'with %d workers: %s'], preferredProfile, maxWorkers, ME1.message);
    end

    retryWorkers = max(1, min(2, maxWorkers));
    if retryWorkers ~= maxWorkers
        try
            poolObj = parpool(preferredProfile, retryWorkers);
            statusMsg = sprintf(['Initial pool startup failed; retry succeeded ' ...
                'with profile=%s, workers=%d.'], preferredProfile, poolObj.NumWorkers);
            return;
        catch ME2
            statusMsg = sprintf('%s\nRetry with fewer workers also failed: %s', ...
                statusMsg, ME2.message);
        end
    end

    poolObj = [];
end
