clear;
close all;

%--------------------------------------------------------------------------
%
% This code forms a simulation of a wideband wireless communications system
% in multipath fading. It simulates: OFDM,  OTFS, coded-OFDM and coded-OTFS
% 
% The following files are required:
% dataGen.m
% multipathChannel.m
% modOFDM.m
% demodOFDM.m
% ISFFT.m
% SFFT.m
% equaliser.m
% plotGraphs.m
%
%--------------------------------------------------------------------------
%
% Author: Bradley Bates
% University of Bristol, UK
% email address: bb16177@bristol.ac.uk
% May 2020
%
% Copyright (c) 2020, Bradley Bates
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Define simulation parameters
%--------------------------------------------------------------------------
M = 16;                         % Modulation alphabet
k = log2(M);                    % Bits/symbol
cpSize = 0.07;                  % OFDM cyclic prefix size
scs = 120e3;                    % Subcarrier spacing Hz (120 kHz for LEO-NTN)
Bw = 100e6;                     % System bandwidth Hz
ofdmSym = 128;                  % No. OFDM symbols / subframe (increased for LEO Doppler)
EbNo = (-3:1:30)';              % Range of energy/bit to noise power ratio
velocity = 26000;               % Velocity of mobile rx relative LEO satellite km/hr
fc = 2e9;                       % Carrier frequency (Hz) S-band for LEO - must match multipathChannel.m
codeRate = 2/4;                 % FEC code rate used
maxIterations = 25;             % Set maximum no. of iterations for LDPC decoder
totalBits = 1e6;                % The approx. total no. of bits simulated
repeats = 1;                    % Number of simulation repetitions

% MP Equalizer Configuration (Raviteja et al., 2018)
mpMaxIter = 10;                 % Max MP iterations (paper: 10 sufficient)
mpDamping = 0.7;                % Damping factor (paper Fig.4: 0.7 optimal)
mpNi = 0;                       % IDI truncation (0=integer Doppler, >0=fractional)


%--------------------------------------------------------------------------
%                    Initialise Simulation Components
%--------------------------------------------------------------------------

% Initialise OFDM Mod/Demod variables
numSC = pow2(ceil(log2(Bw/scs))); % Calc. nearest base-2 no. of OFDM subcarriers
cpLen = floor(cpSize * numSC);    % Calc. cyclic prefix length
numDC = (numSC - 12);             % Calc. number of data carriers

% Initialise the AWGN channel
awgnChannel = comm.AWGNChannel('NoiseMethod','Variance', 'VarianceSource','Input port');
errorRate = comm.ErrorRate('ResetInputPort',true);
errorRate1 = comm.ErrorRate('ResetInputPort',true);

% Initialise the LDPC coder/decoder
parityCheck_matrix = dvbs2ldpc(codeRate);               % Generate DVB-S.2 paraity check matrix
ldpcEncoder = comm.LDPCEncoder(parityCheck_matrix);     % Create encoder system object
ldpcDecoder = comm.LDPCDecoder(parityCheck_matrix);     % Create decoder system object
ldpcDecoder.MaximumIterationCount = maxIterations;      % Assign decoder's maximum iterations
noCodedbits = size(parityCheck_matrix,2);               % Find the Codeword length

% Pre-compute QAM constellation for MP equalizer
qamAlphabet = qammod((0:M-1)', M, 'UnitAveragePower', true);

%--------------------------------------------------------------------------
% LEO Satellite Pilot Pattern Configuration (DD domain)
%--------------------------------------------------------------------------
% Guard bands sized to match LEO channel delay/Doppler spread
% Must match multipathChannel.m parameters (500 ns delay spread cap)
delta_tau_s = 1 / (numSC * scs);                        % Delay resolution (s)
maxDelayspread_s = min(0.5 * cpSize / scs, 500e-9);     % LEO capped delay spread
v_ms = velocity * 1e3 / 3600;                           % Velocity in m/s
fd_hz = v_ms * fc / physconst('LightSpeed');             % Max Doppler shift (Hz)
Ts_sym = (1 + cpSize) / scs;                            % OFDM symbol duration with CP
delta_nu_hz = 1 / (ofdmSym * Ts_sym);                   % Doppler resolution (Hz)

% Doppler pre-compensation from satellite ephemeris (standard in NR-NTN)
% Removes bulk Doppler, leaving only residual scatter spread
fd_precomp_hz = fd_hz;                                      % Pre-comp value (from ephemeris)
residualScatter_hz = 0.05 * fd_hz;                          % Residual scatter after pre-comp
residualScatter_bins = ceil(residualScatter_hz / delta_nu_hz);

maxDelayBins = ceil(maxDelayspread_s / delta_tau_s);            % Physical max delay (bins)

pilotConfig.kp = ceil(ofdmSym/2);                              % Pilot Doppler index (center)
pilotConfig.lp = ceil((Bw/scs - 12)/2);                        % Pilot delay index (center)
pilotConfig.lGuard = 2 * maxDelayBins + 1;                     % Delay guard: 2x for data protection
pilotConfig.kGuard = 2 * residualScatter_bins + 2;             % Doppler guard: residual only (pre-comp)
pilotConfig.maxDelayBins = maxDelayBins;                        % Scan bound: physical delay range
pilotConfig.maxDopplerBins = residualScatter_bins + 1;          % Scan bound: residual Doppler range
pilotConfig.pilotBoostdB = 10;                                  % Pilot power boost (dB)

% Create Vectors for storing error data
berOFDM = zeros(length(EbNo),3); berCOFDM = zeros(length(EbNo),3); berOTFS = zeros(length(EbNo),3); berCOTFS = zeros(length(EbNo),3);
berOTFS_pilot = zeros(length(EbNo),3); berCOTFS_pilot = zeros(length(EbNo),3);
errorStats_coded = zeros(1,3); errorStats_uncoded = zeros(1,3);

for repetition=1:repeats                                % Repeat simulation multiple times with a unqique channel for each repeat
    
    % Generate and Encode data
    [dataIn, dataBits_in, codedData_in, packetSize, numPackets, numCB] = dataGen(k,numDC,ofdmSym,totalBits,codeRate,ldpcEncoder);
    
    % Generate TF-domain channel matrix H(subcarrier, symbol)
    tfGridSize = zeros(numSC, ofdmSym);                              % TF grid dimensions
    [channelTF, chInfo] = multipathChannel(cpSize, scs, tfGridSize, velocity); % TF-domain channel + DD params
    
    % QAM Modulator
    qamTx = qammod(dataIn,M,'InputType','bit','UnitAveragePower',true);    % Apply QAM modulation
    parallelTx = reshape(qamTx,[numDC,ofdmSym*packetSize]);                % Convert to parallel
    
    % Add nulls at index 1
    guardbandTx = [zeros(1,ofdmSym*packetSize); parallelTx];
    % Add 11 nulls around DC
    guardbandTx = [guardbandTx(1:(numDC/2),:); zeros(11,ofdmSym*packetSize); guardbandTx((numDC/2)+1:end,:)];
    
    
%--------------------------------------------------------------------------
%                       OFDM BER Calculation
%--------------------------------------------------------------------------
    
    % Calculate SNR
    snr = EbNo + 10*log10(codeRate*k) + 10*log10(numDC/((numSC)));
    
    % Multicarrier Modulation (OFDM)
    % Apply channel in TF domain, then modulate to time domain
    frameBuffer = guardbandTx;
    txframeBuffer_faded = [];           % Faded time-domain signals
    txframeBuffer_clean = [];           % Clean time-domain signals
    for w = 1:packetSize
        subframeTF = frameBuffer(:,1:ofdmSym);                              % TF-domain subframe
        fadedTF = channelTF .* subframeTF;                                  % Apply channel in TF domain
        ofdmTx_faded = modOFDM(fadedTF, numSC, cpLen, ofdmSym);            % Faded → time domain
        ofdmTx_clean = modOFDM(subframeTF, numSC, cpLen, ofdmSym);         % Clean → time domain
        frameBuffer(:, 1:ofdmSym) = [];
        txframeBuffer_faded = [txframeBuffer_faded; ofdmTx_faded];
        txframeBuffer_clean = [txframeBuffer_clean; ofdmTx_clean];
    end


    % Loop through different values of EbNo
    for m = 1:length(EbNo)
        % Loop through the of packets to be transmitted
        for j = 1:numPackets
            rxframeBuffer = [];                 % Initialise matrix

            % Transmit each subframe individually
            for u = 1:packetSize

                % Extract faded and clean subframes
                fadedSig = txframeBuffer_faded( ((u-1)*numel(ofdmTx_faded)+1) : u*numel(ofdmTx_faded) );
                txSig = txframeBuffer_clean( ((u-1)*numel(ofdmTx_clean)+1) : u*numel(ofdmTx_clean) );

                % AWGN Channel
                release(awgnChannel);
                powerDB = 10*log10(var(fadedSig));            % Calculate Tx signal power
                noiseVar = 10.^(0.1*(powerDB-snr(m)));        % Calculate the noise variance
                rxSig = awgnChannel(fadedSig,noiseVar);       % Pass the signal through a noisy channel

                % Equalisation
                eqSig = equaliser(rxSig,fadedSig,txSig,ofdmSym);

                % Demodulation
                rxSubframe = demodOFDM(eqSig,cpLen,ofdmSym);     % Apply OFDM demodulation
                rxframeBuffer = [rxframeBuffer';rxSubframe']';         % Store demodulated subframe in rx buffer
            end
            % Remove all null carriers
            parallelRx = rxframeBuffer;
            parallelRx((numDC/2)+1:(numDC/2)+11, :) = [];     % Remove nulls around the DC input
            parallelRx(1:1, :) = [];                          % Remove nulls at index 1
            qamRx = reshape(parallelRx,[numel(parallelRx),1]);% Convert to serial
            
            % Uncoded demodulation of entire packet
            dataOut = qamdemod(qamRx,M,'OutputType','bit','UnitAveragePower',true);% Apply QAM demodulation
            codedData_out = randdeintrlv(dataOut,4831);                            % De-interleave data
            codedData_out(numel(codedData_in)+1:end) = [];                         % Remove pad bits
            errorStats_uncoded = errorRate(codedData_in,codedData_out,0);          % Collect error statistics
          
            % Coded demodulation of entire packet
            powerDB = 10*log10(var(qamRx));                                   % Calculate Rx signal power
            noiseVar = 10.^(0.1*(powerDB-(EbNo(m) + 10*log10(codeRate*k) - 10*log10(sqrt(numDC)))));            % Calculate the noise variance
            dataOut = qamdemod(qamRx,M,'OutputType', 'approxllr','UnitAveragePower',true,'NoiseVariance',noiseVar);% Apply QAM demodulation
            codedData_out1 = randdeintrlv(dataOut,4831);                      % De-interleave data
            codedData_out1(numel(codedData_in)+1:end) = [];                   % Remove pad bits
            
            % Decode individual code blocks
            dataBits_out = [];                                                % Initialise matrix
            dataOut_buffer = codedData_out1;
            for q = 1:numCB
                dataBits_out = [dataBits_out;ldpcDecoder(dataOut_buffer(1:noCodedbits))]; % Decode data & add it to the data bits out matrix
                dataOut_buffer(1:noCodedbits) = [];                                       % Delete decoded data from buffer
            end
            dataBits_out = double(dataBits_out);                              % Convert to a double compatible w/ errorStats
            errorStats_coded = errorRate1(dataBits_in,dataBits_out,0);        % Collect error statistics
            
        end
        berOFDM(m,:) = errorStats_uncoded;                                  % Save uncoded BER data
        berCOFDM(m,:) = errorStats_coded;                                   % Save  coded BER data
        errorStats_uncoded = errorRate(codedData_in,codedData_out,1);       % Reset the error rate calculator
        errorStats_coded = errorRate1(dataBits_in,dataBits_out,1);          % Reset the error rate calculator
        
    end
    
    
%--------------------------------------------------------------------------
%                       OTFS BER Calculation
%--------------------------------------------------------------------------
    
    % Calculate SNR
    snr = EbNo + 10*log10(codeRate*k) + 10*log10(numDC/((numSC))) + 10*log10(sqrt(ofdmSym));
    
    % Multicarrier Modulation (OTFS)
    % Apply channel in TF domain between ISFFT and OFDM mod
    frameBuffer = guardbandTx;
    txframeBuffer_faded = [];
    txframeBuffer_clean = [];
    for w = 1:packetSize
        otfsTF = ISFFT(frameBuffer(:,1:ofdmSym));                  % DD → TF domain
        fadedTF = channelTF .* otfsTF;                              % Apply channel in TF domain
        ofdmTx_faded = modOFDM(fadedTF, numSC, cpLen, ofdmSym);    % Faded → time domain
        ofdmTx_clean = modOFDM(otfsTF, numSC, cpLen, ofdmSym);     % Clean → time domain
        frameBuffer(:, 1:ofdmSym) = [];
        txframeBuffer_faded = [txframeBuffer_faded; ofdmTx_faded];
        txframeBuffer_clean = [txframeBuffer_clean; ofdmTx_clean];
    end

    % Loop through different values of EbNo
    for m = 1:length(EbNo)
        % Loop through the of packets to be transmitted
        for j = 1:numPackets
            rxframeBuffer = [];                 % Initialise matrix

            % Transmit each subframe individually
            for u = 1:packetSize

                % Extract faded and clean subframes
                fadedSig = txframeBuffer_faded( ((u-1)*numel(ofdmTx_faded)+1) : u*numel(ofdmTx_faded) );
                txSig = txframeBuffer_clean( ((u-1)*numel(ofdmTx_clean)+1) : u*numel(ofdmTx_clean) );

                % AWGN Channel
                release(awgnChannel);
                powerDB = 10*log10(var(fadedSig));            % Calculate Tx signal power
                noiseVar = 10.^(0.1*(powerDB-snr(m)));        % Calculate the noise variance
                rxSig = awgnChannel(fadedSig,noiseVar);       % Pass the signal through a noisy channel

                % Equalisation
                eqSig = equaliser(rxSig,fadedSig,txSig,ofdmSym);

                % Demodulation
                otfsRx = demodOFDM(eqSig,cpLen,ofdmSym);     % Apply OFDM demodulation
                rxSubframe = SFFT(otfsRx);                      % Apply OTFS demodulation
                rxframeBuffer = [rxframeBuffer';rxSubframe']';     % Store demodulated subframe in rx buffer
            end
            % Remove all null carriers
            parallelRx = rxframeBuffer;
            parallelRx((numDC/2)+1:(numDC/2)+11, :) = [];         % Remove nulls around the DC input
            parallelRx(1:1, :) = [];                              % Remove nulls at index 1
            qamRx = reshape(parallelRx,[numel(parallelRx),1]);    % Convert to serial
            
            % Uncoded demodulation of entire packet
            dataOut = qamdemod(qamRx,M,'OutputType','bit','UnitAveragePower',true);% Apply QAM demodulation
            codedData_out = randdeintrlv(dataOut,4831);                            % De-interleave data
            codedData_out(numel(codedData_in)+1:end) = [];                         % Remove pad bits
            errorStats_uncoded = errorRate(codedData_in,codedData_out,0);          % Collect error statistics
            

            % Coded demodulation of entire packet
            powerDB = 10*log10(var(qamRx));                                   % Calculate Rx signal power
            noiseVar = 10.^(0.1*(powerDB-(EbNo(m) + 10*log10(codeRate*k) - 10*log10(sqrt(numDC)))));            % Calculate the noise variance
            dataOut = qamdemod(qamRx,M,'OutputType', 'approxllr','UnitAveragePower',true,'NoiseVariance',noiseVar);% Apply QAM demodulation
            codedData_out1 = randdeintrlv(dataOut,4831);                      % De-interleave data
            codedData_out1(numel(codedData_in)+1:end) = [];                   % Remove pad bits
            
            % Decode individual code blocks
            dataBits_out = [];                                                % Initialise matrix
            dataOut_buffer = codedData_out1;
            for q = 1:numCB
                dataBits_out = [dataBits_out;ldpcDecoder(dataOut_buffer(1:noCodedbits))]; % Decode data & add it to the data bits out matrix
                dataOut_buffer(1:noCodedbits) = [];                                       % Delete decoded data from buffer
            end
            dataBits_out = double(dataBits_out);                              % Convert to a double compatible w/ errorStats
            errorStats_coded = errorRate1(dataBits_in,dataBits_out,0);     % Collect error statistics
            
        end
        berOTFS(m,:) = errorStats_uncoded;                                  % Save uncoded BER data
        berCOTFS(m,:) = errorStats_coded;                                   % Save coded BER data
        errorStats_uncoded = errorRate(codedData_in,codedData_out,1);       % Reset the error rate calculator
        errorStats_coded = errorRate1(dataBits_in,dataBits_out,1);          % Reset the error rate calculator

    end


%--------------------------------------------------------------------------
%              OTFS with DD-Domain Pilot Pattern (LEO Satellite)
%--------------------------------------------------------------------------

    % Calculate SNR (same as OTFS)
    snr = EbNo + 10*log10(codeRate*k) + 10*log10(numDC/((numSC))) + 10*log10(sqrt(ofdmSym));

    % --- Transmitter: Place pilot + data on DD grid ---
    % For each subframe, create DD grid with pilot pattern
    frameBuffer_coded = guardbandTx;    % Use same guard-banded QAM data

    % Store pilot info for all subframes (same pattern each subframe)
    pilotConfig.lp = ceil(size(guardbandTx,1)/2);  % Update to match actual grid rows
    pilotConfig.kp = ceil(ofdmSym/2);

    txDDGrids = cell(packetSize, 1);
    rxDDGrids_faded = cell(packetSize, 1);
    ddGridSize = [size(guardbandTx,1), ofdmSym];

    for w = 1:packetSize
        subframe = frameBuffer_coded(:, 1:ofdmSym);
        frameBuffer_coded(:, 1:ofdmSym) = [];

        % Extract data symbols from subframe (column-major order)
        dataSyms = subframe(:);

        % Create DD grid with pilot pattern
        [ddGrid, dataIdx, pilotIdx, guardIdx, pilotInfoTx] = ...
            pilotPatternDD(dataSyms, size(subframe,1), ofdmSym, pilotConfig);

        % Pass scan bounds to pilotInfoTx for channel estimation
        pilotInfoTx.maxDelayBins = pilotConfig.maxDelayBins;
        pilotInfoTx.maxDopplerBins = pilotConfig.maxDopplerBins;

        % Apply channel in DD domain with Doppler pre-compensation
        txDDGrids{w} = ddGrid;
        rxDDGrids_faded{w} = applyChannelDD(ddGrid, chInfo, scs, cpSize, fd_precomp_hz);
    end

    % --- Receiver: Add noise in DD domain and estimate/equalize ---
    for m = 1:length(EbNo)
        for j = 1:numPackets
            rxframeBuffer = [];

            for u = 1:packetSize
                rxDD_faded = rxDDGrids_faded{u};

                % Add AWGN directly in DD domain
                % (SFFT/ISFFT are unitary transforms, so AWGN statistics are preserved)
                sigPower = 10*log10(mean(abs(rxDD_faded(:)).^2));
                noiseVar = 10^(0.1*(sigPower - snr(m)));
                noise = sqrt(noiseVar/2) * (randn(size(rxDD_faded)) + 1j*randn(size(rxDD_faded)));
                rxDD = rxDD_faded + noise;

                % DD-domain channel estimation using pilot
                [hEst, delayEst, dopplerEst, navInfo] = ddChannelEstimate(rxDD, pilotInfoTx);

                % --- Pilot Interference Cancellation (PIC) ---
                % Reconstruct and subtract pilot response from received grid
                % navInfo.pathGains = pilotAmp * h_i (includes pilot scaling)
                pilotAmp = txDDGrids{u}(pilotInfoTx.lp, pilotInfoTx.kp);
                pilotResponse = zeros(ddGridSize);
                for pp = 1:navInfo.numPathsDetected
                    lr = pilotInfoTx.lp + navInfo.pathDelays(pp);
                    kc = pilotInfoTx.kp + navInfo.pathDopplers(pp);
                    if lr >= 1 && lr <= ddGridSize(1) && kc >= 1 && kc <= ddGridSize(2)
                        pilotResponse(lr, kc) = navInfo.pathGains(pp);
                    end
                end
                rxDD_clean = rxDD - pilotResponse;

                % --- Build sparse DD channel matrix from estimated paths ---
                % Normalize gains: remove pilot amplitude to get pure channel h_i
                channelGains = navInfo.pathGains / pilotAmp;
                H_dd = buildDDChannelMatrix(navInfo.pathDelays, navInfo.pathDopplers, ...
                    channelGains, ddGridSize(1), ddGridSize(2), mpNi);

                % --- Message Passing equalization (Raviteja 2018) ---
                mpCfg.maxIter = mpMaxIter;
                mpCfg.damping = mpDamping;
                [eqDataSyms, ~] = otfsEqualizerMP(rxDD_clean(:), H_dd, dataIdx, ...
                    qamAlphabet, noiseVar, mpCfg);

                % Reconstruct subframe from equalized data at data positions
                rxSubframe = zeros(ddGridSize);
                if length(eqDataSyms) >= length(dataIdx)
                    rxSubframe(dataIdx) = eqDataSyms(1:length(dataIdx));
                end

                rxframeBuffer = [rxframeBuffer'; rxSubframe']';
            end

            % Remove guard band nulls (same as original OTFS path)
            parallelRx = rxframeBuffer;
            parallelRx((numDC/2)+1:(numDC/2)+11, :) = [];
            parallelRx(1:1, :) = [];
            qamRx = reshape(parallelRx, [numel(parallelRx), 1]);

            % Uncoded demodulation
            dataOut = qamdemod(qamRx, M, 'OutputType', 'bit', 'UnitAveragePower', true);
            codedData_out = randdeintrlv(dataOut, 4831);
            codedData_out(numel(codedData_in)+1:end) = [];
            errorStats_uncoded = errorRate(codedData_in, codedData_out, 0);

            % Coded demodulation
            powerDB = 10*log10(var(qamRx));
            noiseVar = 10.^(0.1*(powerDB-(EbNo(m) + 10*log10(codeRate*k) - 10*log10(sqrt(numDC)))));
            dataOut = qamdemod(qamRx, M, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', noiseVar);
            codedData_out1 = randdeintrlv(dataOut, 4831);
            codedData_out1(numel(codedData_in)+1:end) = [];

            dataBits_out = [];
            dataOut_buffer = codedData_out1;
            for q = 1:numCB
                dataBits_out = [dataBits_out; ldpcDecoder(dataOut_buffer(1:noCodedbits))];
                dataOut_buffer(1:noCodedbits) = [];
            end
            dataBits_out = double(dataBits_out);
            errorStats_coded = errorRate1(dataBits_in, dataBits_out, 0);
        end

        berOTFS_pilot(m,:) = errorStats_uncoded;
        berCOTFS_pilot(m,:) = errorStats_coded;
        errorStats_uncoded = errorRate(codedData_in, codedData_out, 1);
        errorStats_coded = errorRate1(dataBits_in, dataBits_out, 1);
    end

    % Print navigation info from last subframe
    fprintf('\n--- LEO Navigation Info (from DD-domain pilot) ---\n');
    fprintf('Detected %d channel paths\n', navInfo.numPathsDetected);
    fprintf('Pilot SNR: %.1f dB\n', navInfo.snrPilot);
    fprintf('Pilot overhead: %.1f%%\n', pilotInfoTx.overheadPercent);
    fprintf('Effective pilot boost: %.1f dB\n', pilotInfoTx.effectivePilotBoostdB);
    for pp = 1:length(navInfo.pathDelays)
        fprintf('  Path %d: delay=%+d bins, Doppler=%+d bins, |gain|=%.3f\n', ...
            pp, navInfo.pathDelays(pp), navInfo.pathDopplers(pp), abs(navInfo.pathGains(pp)));
    end
    fprintf('---------------------------------------------------\n');

    %----------------------------------------------------------------------
    %   Velocity & Range Estimation from DD-Domain Pilot Detection
    %----------------------------------------------------------------------
    % Physical resolutions of the OTFS DD grid
    c = physconst('LightSpeed');            % Speed of light (m/s)
    Ts_sym = (1 + cpSize) / scs;            % OFDM symbol duration with CP (s)
    N_delay = size(rxDD, 1);                % Number of delay bins (= numSC)
    M_doppler = size(rxDD, 2);              % Number of Doppler bins (= ofdmSym)

    delta_tau = 1 / (N_delay * scs);        % Delay resolution (s/bin)
    delta_nu  = 1 / (M_doppler * Ts_sym);   % Doppler resolution (Hz/bin)

    fprintf('\n=== Velocity & Range Estimation ===\n');
    fprintf('System parameters:\n');
    fprintf('  Carrier frequency fc     = %.2f GHz\n', fc/1e9);
    fprintf('  Subcarrier spacing df    = %.0f kHz\n', scs/1e3);
    fprintf('  Delay bins N             = %d\n', N_delay);
    fprintf('  Doppler bins M           = %d\n', M_doppler);
    fprintf('  Delay resolution  d_tau  = %.2f ns  (%.2f m)\n', delta_tau*1e9, delta_tau*c);
    fprintf('  Doppler resolution d_nu  = %.2f Hz  (%.2f m/s)\n', delta_nu, delta_nu*c/fc);

    % --- Identify LoS path (strongest gain) ---
    [~, losIdx] = max(abs(navInfo.pathGains));
    los_delay_bin   = navInfo.pathDelays(losIdx);
    los_doppler_bin = navInfo.pathDopplers(losIdx);

    % Use fractional Doppler/delay from parabolic interpolation
    los_frac_doppler = navInfo.losFracDoppler;  % Fractional Doppler (bins)
    los_frac_delay   = navInfo.losFracDelay;    % Fractional delay (bins)

    % Convert LoS path to physical values (using fractional estimates)
    % Add back the pre-compensated bulk Doppler (from ephemeris)
    % NOTE: SFFT convention maps positive Doppler to negative DD bins
    %       Physical Doppler = -DD_Doppler_shift * delta_nu
    tau_los = los_frac_delay * delta_tau;                   % Propagation delay (s)
    nu_residual = -los_frac_doppler * delta_nu;             % Residual Doppler (Hz) — negate for SFFT sign
    nu_los  = nu_residual + fd_precomp_hz;                  % Total Doppler = residual + bulk (Hz)
    range_los = tau_los * c;                                % One-way range (m)
    v_est = nu_los * c / fc;                                % Estimated radial velocity (m/s)
    v_est_kmh = v_est * 3.6;                                % Convert to km/hr

    % True velocity for comparison
    v_true_kmh = velocity;                      % From simulation parameter (km/hr)
    v_true_ms  = velocity * 1e3 / 3600;         % m/s
    fd_true = v_true_ms * fc / c;               % True max Doppler shift (Hz)

    fprintf('\n  Doppler pre-compensation = %.2f Hz (from ephemeris)\n', fd_precomp_hz);
    fprintf('  Residual Doppler guard   = %d bins (kGuard)\n', pilotConfig.kGuard);
    fprintf('\n--- LoS Path (strongest, Path %d) ---\n', losIdx);
    fprintf('  Delay   = %+d bins (frac: %+.2f)  =>  tau = %.2f ns  =>  range = %.2f m\n', ...
        los_delay_bin, los_frac_delay, tau_los*1e9, range_los);
    fprintf('  Doppler = %+d bins (frac: %+.2f)  =>  nu_residual = %.2f Hz\n', ...
        los_doppler_bin, los_frac_doppler, nu_residual);
    fprintf('  Total Doppler = residual + precomp = %.2f + %.2f = %.2f Hz\n', ...
        nu_residual, fd_precomp_hz, nu_los);
    fprintf('  v_est = %.2f m/s (%.2f km/hr)\n', v_est, v_est_kmh);
    fprintf('\n--- Comparison with True Velocity ---\n');
    fprintf('  True velocity       = %.2f m/s  (%.2f km/hr)\n', v_true_ms, v_true_kmh);
    fprintf('  True max Doppler fd = %.2f Hz\n', fd_true);
    fprintf('  Estimated velocity  = %.2f m/s  (%.2f km/hr)\n', v_est, v_est_kmh);
    fprintf('  Velocity error      = %.2f m/s  (%.2f km/hr)\n', abs(v_est - v_true_ms), abs(v_est_kmh - v_true_kmh));
    fprintf('  Relative error      = %.1f%%\n', abs(v_est - v_true_ms)/v_true_ms * 100);

    % --- All paths: physical delay and Doppler ---
    fprintf('\n--- All Detected Paths (physical values) ---\n');
    fprintf('  %-6s  %-12s  %-12s  %-14s  %-14s  %-10s\n', ...
        'Path', 'Delay(ns)', 'Range(m)', 'Doppler(Hz)', 'Velocity(m/s)', '|Gain|');
    for pp = 1:length(navInfo.pathDelays)
        tau_pp = navInfo.pathDelays(pp) * delta_tau;
        nu_pp  = -navInfo.pathDopplers(pp) * delta_nu + fd_precomp_hz;  % Negate for SFFT sign + add bulk Doppler
        range_pp = tau_pp * c;
        v_pp = nu_pp * c / fc;
        fprintf('  %-6d  %-12.2f  %-12.2f  %-14.2f  %-14.2f  %-10.3f\n', ...
            pp, tau_pp*1e9, range_pp, nu_pp, v_pp, abs(navInfo.pathGains(pp)));
    end
    fprintf('==========================================\n');

end

%--------------------------------------------------------------------------
%                           Figures
%-------------------------------------------------------------------------- 

% Plot BER / EbNo curves (including OTFS-Pilot results)
plotGraphs(berOFDM, berCOFDM, berOTFS, berCOTFS, M, numSC, EbNo, berOTFS_pilot, berCOTFS_pilot);

% Plot DD-domain pilot pattern visualization
figure;
ddVis = zeros(size(guardbandTx,1), ofdmSym);
ddVis(dataIdx) = 0.3;                  % Data positions (gray)
ddVis(guardIdx) = 0;                   % Guard band (black)
ddVis(pilotIdx) = 1;                   % Pilot (bright)
imagesc(ddVis);
colormap(hot);
colorbar;
title('Delay-Doppler Grid: Pilot Pattern with Guard Bands');
xlabel('Doppler Index (OFDM symbols)');
ylabel('Delay Index (Subcarriers)');
hold on;
% Draw guard band rectangle
rectangle('Position', [pilotConfig.kp-pilotConfig.kGuard-0.5, ...
    pilotConfig.lp-pilotConfig.lGuard-0.5, ...
    2*pilotConfig.kGuard+1, 2*pilotConfig.lGuard+1], ...
    'EdgeColor', 'c', 'LineWidth', 2, 'LineStyle', '--');
legend('Guard Band Boundary');
hold off;

