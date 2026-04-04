function [hEst, delayEst, dopplerEst, navInfo] = ddChannelEstimate(rxGrid, pilotInfo, fracCancel)

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
%   fracCancel mode (P0-1 enhancement):
%   When fracCancel=true, uses Quinn fractional estimates + Dirichlet
%   kernel to reconstruct the full 2D spreading response of each detected
%   path before pilot cancellation. This closes the model-algorithm gap
%   between the fractional DD I/O model and the estimation/cancellation.
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
% fracCancel        (optional) Enable fractional-aware pilot cancellation.
%                   Default = false (legacy delta cancellation).
%
%--------------------------------------------------------------------------
% Function returns:
%
% hEst              N x M estimated DD domain channel matrix
% delayEst          Estimated delay indices of channel paths
% dopplerEst        Estimated Doppler indices of channel paths
% navInfo           Struct with navigation-relevant measurements
%                     .pilotResponse  - (when fracCancel=true) full N x M
%                       reconstructed pilot response including Dirichlet
%                       spreading, ready for subtraction from rxGrid
%
%--------------------------------------------------------------------------

if nargin < 3, fracCancel = false; end

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
scanRows = lp : min(N, lp + maxDelayBins);
scanCols = max(1, kp - maxDopplerBins) : min(M, kp + maxDopplerBins);

%% Detection threshold
scanRegion = abs(rxGrid(scanRows, scanCols));
pilotAmp = max(scanRegion(:));
noiseFloor = median(scanRegion(:));
threshRelativedB = -25;
threshPilot = pilotAmp * 10^(threshRelativedB/20);
threshNoise = 3.0 * noiseFloor;
threshold = max(threshPilot, threshNoise);

%% Build DD channel impulse response
hEst = zeros(N, M);

maxPaths = 10;
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

%% Per-path fractional estimation via Quinn interpolator
% For each detected path, estimate fractional delay and Doppler offsets
numDetected = length(delayEst);
fracDelays = zeros(numDetected, 1);    % fractional delay for each path
fracDopplers = zeros(numDetected, 1);  % fractional Doppler for each path

for pp = 1:numDetected
    row_pp = lp + delayEst(pp);
    col_pp = kp + dopplerEst(pp);

    % Quinn along Doppler (column) dimension
    delta_k = 0;
    if col_pp > 1 && col_pp < M
        Xm1 = rxGrid(row_pp, col_pp - 1);
        X0  = rxGrid(row_pp, col_pp);
        Xp1 = rxGrid(row_pp, col_pp + 1);
        if abs(X0) > 1e-15
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
        end
    end
    fracDopplers(pp) = delta_k;

    % Quinn along delay (row) dimension
    delta_l = 0;
    if row_pp > 1 && row_pp < N
        Xm1 = rxGrid(row_pp - 1, col_pp);
        X0  = rxGrid(row_pp, col_pp);
        Xp1 = rxGrid(row_pp + 1, col_pp);
        if abs(X0) > 1e-15
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
        end
    end
    fracDelays(pp) = delta_l;
end

%% LoS path selection (for navigation)
maxAmp = max(abs(pathGains));
strongMask = abs(pathGains) > maxAmp * 10^(-6/20);  % Within 6 dB
strongIdx = find(strongMask);
[~, minPos] = min(abs(delayEst(strongIdx)));
losIdx = strongIdx(minPos);

fracDopplerShift = dopplerEst(losIdx) + fracDopplers(losIdx);
fracDelayShift = delayEst(losIdx) + fracDelays(losIdx);

%% Fractional-aware pilot response reconstruction
% Reconstruct the full pilot spreading response using Dirichlet kernels
% for each detected path's fractional delay and Doppler offsets.
if fracCancel
    Ns = 2;  % Dirichlet spreading half-width (must match applyChannelDD)
    pilotResponse = zeros(N, M);
    l_idx = (0:N-1)';

    for pp = 1:numDetected
        hi_pp = pathGains(pp);  % This is already h_i * x_p from the guard region

        li_int = delayEst(pp);
        ki_int = dopplerEst(pp);
        lambda = fracDelays(pp);   % fractional part [-0.5, 0.5]
        kappa = fracDopplers(pp);  % fractional part [-0.5, 0.5]

        % Phase correction
        ki_full = ki_int + kappa;
        phase_col = exp(1j * 2 * pi * l_idx * ki_full / (N * M));

        nearIntDelay = (abs(lambda) < 1e-6);
        nearIntDoppler = (abs(kappa) < 1e-6);

        if nearIntDelay && nearIntDoppler
            % No spreading — place as delta
            r_idx = mod(lp + li_int - 1, N) + 1;
            c_idx = mod(kp + ki_int - 1, M) + 1;
            pilotResponse(r_idx, c_idx) = pilotResponse(r_idx, c_idx) + hi_pp;
        else
            % Reconstruct Dirichlet spreading
            for p = -Ns:Ns
                if abs(lambda - p) < 1e-10
                    alpha_p = 1;
                else
                    alpha_p = (1/N) * (1 - exp(1j*2*pi*(lambda - p))) / ...
                              (1 - exp(1j*2*pi*(lambda - p)/N));
                end
                if abs(alpha_p) < 1e-15, continue; end

                for q = -Ns:Ns
                    if abs(kappa - q) < 1e-10
                        beta_q = 1;
                    else
                        beta_q = (1/M) * (1 - exp(1j*2*pi*(kappa - q))) / ...
                                 (1 - exp(1j*2*pi*(kappa - q)/M));
                    end
                    if abs(beta_q) < 1e-15, continue; end

                    r_idx = mod(lp + li_int + p - 1, N) + 1;
                    c_idx = mod(kp + ki_int + q - 1, M) + 1;
                    pilotResponse(r_idx, c_idx) = pilotResponse(r_idx, c_idx) + ...
                        hi_pp * alpha_p * beta_q * phase_col(r_idx);
                end
            end
        end
    end
else
    % Legacy delta cancellation: simple impulse at detected bins
    pilotResponse = zeros(N, M);
    for pp = 1:numDetected
        lr = mod(lp + delayEst(pp) - 1, N) + 1;
        kc = mod(kp + dopplerEst(pp) - 1, M) + 1;
        pilotResponse(lr, kc) = pathGains(pp);
    end
end

%% Navigation information extraction
pilotEnergy = pilotAmp^2;
noiseEst = median(abs(rxGrid(:)).^2) * 0.5;

navInfo.pathDelays = delayEst;
navInfo.pathDopplers = dopplerEst;
navInfo.pathGains = pathGains;
navInfo.fracDelays = fracDelays;
navInfo.fracDopplers = fracDopplers;
navInfo.pilotPower = pilotEnergy;
navInfo.noiseEst = noiseEst;
if noiseEst > 0
    navInfo.snrPilot = 10*log10(pilotEnergy / noiseEst);
else
    navInfo.snrPilot = Inf;
end
navInfo.numPathsDetected = numDetected;
navInfo.losFracDoppler = fracDopplerShift;
navInfo.losFracDelay = fracDelayShift;
navInfo.pilotResponse = pilotResponse;

end
