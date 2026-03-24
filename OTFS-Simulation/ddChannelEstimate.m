function [hEst, delayEst, dopplerEst, navInfo] = ddChannelEstimate(rxGrid, pilotInfo)

%--------------------------------------------------------------------------
%
%   Estimates the channel in the Delay-Doppler domain using the received
%   pilot response. Also extracts navigation parameters (delay/Doppler)
%   for LEO satellite integrated communication and navigation.
%
%   Method: The transmitted pilot is a single impulse at (lp, kp).
%   After passing through an L-path channel, the received DD grid shows
%   the channel response as shifted/scaled copies of the pilot.
%   The channel taps h_i at (l_i, k_i) can be read directly from the
%   guard region of the received grid.
%
%   Scanning is restricted to physically valid regions:
%   - Delay: [0, maxDelayBins] (channel delays are non-negative)
%   - Doppler: [-maxDopplerBins, +maxDopplerBins] (residual after pre-comp)
%   This prevents data contamination from being detected as false paths.
%
%--------------------------------------------------------------------------
% Input arguments:
%
% rxGrid            N x M received DD grid (after SFFT demodulation)
% pilotInfo         Struct from pilotPatternDD with pilot position info:
%                     .lp, .kp       - pilot position
%                     .lGuard, .kGuard - guard band sizes (for data placement)
%                     .maxDelayBins  - max physical delay (bins), scan bound
%                     .maxDopplerBins - max residual Doppler (bins), scan bound
%
%--------------------------------------------------------------------------
% Function returns:
%
% hEst              N x M estimated DD domain channel matrix
% delayEst          Estimated delay indices of channel paths
% dopplerEst        Estimated Doppler indices of channel paths
% navInfo           Struct with navigation-relevant measurements
%
%--------------------------------------------------------------------------

[N, M] = size(rxGrid);
kp = pilotInfo.kp;
lp = pilotInfo.lp;
kGuard = pilotInfo.kGuard;
lGuard = pilotInfo.lGuard;

% Physical scan bounds (restrict to valid channel response region)
if isfield(pilotInfo, 'maxDelayBins')
    maxDelayBins = pilotInfo.maxDelayBins;
else
    maxDelayBins = lGuard;   % Fallback: scan full guard
end
if isfield(pilotInfo, 'maxDopplerBins')
    maxDopplerBins = pilotInfo.maxDopplerBins;
else
    maxDopplerBins = kGuard;  % Fallback: scan full guard
end

%% Define scan region (physically valid channel response area only)
% Delay: channel delays are non-negative, so scan [lp, lp + maxDelayBins]
% Doppler: after pre-comp, residual is bounded by ±maxDopplerBins
scanRows = lp : min(N, lp + maxDelayBins);
scanCols = max(1, kp - maxDopplerBins) : min(M, kp + maxDopplerBins);

%% Detection threshold: pilot-relative
% Real channel responses have amplitude = pilot_amplitude * channel_gain.
% Use a threshold relative to the pilot to reject data contamination.
% With pilot boost of ~10 dB, pilot amplitude ≈ 3.16x data amplitude.
% Threshold at -8 dB below pilot catches paths with |h| > 0.4 (rel to LoS)
% while rejecting data contamination (amplitude ~ 1/pilot_boost).
%
% Use max amplitude in scan region (not fixed (lp,kp)) because bulk
% propagation delay can shift the pilot response away from the nominal
% position, leaving only noise at (lp,kp).
scanRegion = abs(rxGrid(scanRows, scanCols));
pilotAmp = max(scanRegion(:));
% Adaptive threshold: use noise floor estimate to detect weak scatter paths
% Noise floor from median of scan region (robust against signal peaks)
noiseFloor = median(scanRegion(:));
% Threshold = max(relative to pilot, 3x noise floor)
% This detects paths down to -25 dB below pilot while avoiding noise
threshRelativedB = -25;
threshPilot = pilotAmp * 10^(threshRelativedB/20);
threshNoise = 3.0 * noiseFloor;
threshold = max(threshPilot, threshNoise);

%% Build DD channel impulse response
hEst = zeros(N, M);

% Scan the valid region for significant taps
% Collect all candidates above threshold, then keep strongest
maxPaths = 10;   % Limit to prevent noise over-detection
candidateDelays = [];
candidateDopplers = [];
candidateGains = [];
candidateAmps = [];

for r = 1:length(scanRows)
    for c = 1:length(scanCols)
        lr = scanRows(r);
        kc = scanCols(c);
        val = rxGrid(lr, kc);

        if abs(val) > threshold
            candidateDelays = [candidateDelays; lr - lp];
            candidateDopplers = [candidateDopplers; kc - kp];
            candidateGains = [candidateGains; val];
            candidateAmps = [candidateAmps; abs(val)];
        end
    end
end

% Keep only the strongest paths (up to maxPaths)
if length(candidateAmps) > maxPaths
    [~, sortIdx] = sort(candidateAmps, 'descend');
    keepIdx = sortIdx(1:maxPaths);
    candidateDelays = candidateDelays(keepIdx);
    candidateDopplers = candidateDopplers(keepIdx);
    candidateGains = candidateGains(keepIdx);
end

delayEst = candidateDelays;
dopplerEst = candidateDopplers;
pathGains = candidateGains;

% Place in DD channel matrix
for idx = 1:length(delayEst)
    delayIdx = mod(delayEst(idx), N) + 1;
    dopplerIdx = mod(dopplerEst(idx), M) + 1;
    hEst(delayIdx, dopplerIdx) = pathGains(idx);
end

%% If no paths detected, use the pilot position as single-tap channel
if isempty(delayEst)
    hEst(1, 1) = rxGrid(lp, kp);
    delayEst = 0;
    dopplerEst = 0;
    pathGains = rxGrid(lp, kp);
end

%% Fractional Doppler estimation for LoS path (Quinn estimator)
% Quinn's second estimator uses complex DFT bins for optimal accuracy.
% LoS selection: among strong paths (within 6 dB of strongest), pick
% the one with smallest delay. This avoids both:
%   1. Dirichlet sidelobes (weaker, filtered out by 6 dB threshold)
%   2. Scatter paths during shadow fading (larger delay than LoS)
maxAmp = max(abs(pathGains));
strongMask = abs(pathGains) > maxAmp * 10^(-6/20);  % Within 6 dB
strongIdx = find(strongMask);
[~, minPos] = min(abs(delayEst(strongIdx)));
losIdx = strongIdx(minPos);
losRow = lp + delayEst(losIdx);   % row index in rxGrid
losCol = kp + dopplerEst(losIdx); % col index in rxGrid

% Quinn estimator along Doppler (column) dimension
fracDopplerShift = dopplerEst(losIdx);  % default: integer bin
if losCol > 1 && losCol < M
    Xm1 = rxGrid(losRow, losCol - 1);
    X0  = rxGrid(losRow, losCol);
    Xp1 = rxGrid(losRow, losCol + 1);

    ap = real(Xp1 / X0);
    am = real(Xm1 / X0);
    dp = -ap / (1 - ap);
    dm = am / (1 - am);

    if abs(dp) < abs(dm)
        delta_k = dp;
    else
        delta_k = dm;
    end

    delta_k = max(-0.5, min(0.5, delta_k));
    fracDopplerShift = dopplerEst(losIdx) + delta_k;
end

% Quinn estimator along delay (row) dimension
fracDelayShift = delayEst(losIdx);  % default: integer bin
if losRow > 1 && losRow < N
    Xm1 = rxGrid(losRow - 1, losCol);
    X0  = rxGrid(losRow, losCol);
    Xp1 = rxGrid(losRow + 1, losCol);

    ap = real(Xp1 / X0);
    am = real(Xm1 / X0);
    dp = -ap / (1 - ap);
    dm = am / (1 - am);

    if abs(dp) < abs(dm)
        delta_l = dp;
    else
        delta_l = dm;
    end

    delta_l = max(-0.5, min(0.5, delta_l));
    fracDelayShift = delayEst(losIdx) + delta_l;
end

%% Navigation information extraction
pilotEnergy = pilotAmp^2;  % Use actual peak, not fixed (lp,kp)
noiseEst = median(abs(rxGrid(:)).^2) * 0.5;  % Noise from full grid

navInfo.pathDelays = delayEst;
navInfo.pathDopplers = dopplerEst;
navInfo.pathGains = pathGains;
navInfo.pilotPower = pilotEnergy;
navInfo.noiseEst = noiseEst;
if noiseEst > 0
    navInfo.snrPilot = 10*log10(pilotEnergy / noiseEst);
else
    navInfo.snrPilot = Inf;
end
navInfo.numPathsDetected = length(delayEst);
navInfo.losFracDoppler = fracDopplerShift;
navInfo.losFracDelay = fracDelayShift;

end
