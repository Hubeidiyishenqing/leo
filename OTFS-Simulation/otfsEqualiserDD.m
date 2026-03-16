function [eqGrid, dataSymbols] = otfsEqualiserDD(rxGrid, hEst, dataIdx, pilotInfo)

%--------------------------------------------------------------------------
%
%   Equalizes the received OTFS signal in the Delay-Doppler domain
%   using the estimated DD channel response.
%
%   Method: Constructs the effective DD channel matrix and performs
%   MMSE or ZF equalization to recover data symbols.
%
%--------------------------------------------------------------------------
% Input arguments:
%
% rxGrid            N x M received DD grid (after SFFT demodulation)
% hEst              N x M estimated DD channel matrix (from ddChannelEstimate)
% dataIdx           Linear indices of data positions (from pilotPatternDD)
% pilotInfo         Struct with pilot metadata
%
%--------------------------------------------------------------------------
% Function returns:
%
% eqGrid            N x M equalized DD grid
% dataSymbols       Column vector of equalized data symbols (at dataIdx)
%
%--------------------------------------------------------------------------

[N, M] = size(rxGrid);

%% Method: Frequency-domain MMSE equalization
% Convert DD channel to TF domain, equalize there, convert back.
% This is more computationally efficient than building the full
% N*M x N*M DD domain channel matrix.

% Transform received DD grid to TF domain
rxTF = sqrt(N/M) * fft(ifft(rxGrid, [], 1), [], 2);   % ISFFT

% Transform DD channel to TF domain
hTF = sqrt(N/M) * fft(ifft(hEst, [], 1), [], 2);      % ISFFT

% Estimate noise variance from guard region
kp = pilotInfo.kp;
lp = pilotInfo.lp;
kGuard = pilotInfo.kGuard;
lGuard = pilotInfo.lGuard;
guardRows = max(1, lp - lGuard) : min(N, lp + lGuard);
guardCols = max(1, kp - kGuard) : min(M, kp + kGuard);

% Noise from guard positions that are not near detected paths
guardVals = rxGrid(guardRows, guardCols);
% Use the weakest values in guard region as noise estimate
sortedGuard = sort(abs(guardVals(:)).^2);
numNoisesamples = max(1, floor(length(sortedGuard) * 0.3));
noiseVar = mean(sortedGuard(1:numNoisesamples));
if noiseVar == 0
    noiseVar = 1e-10;  % Prevent division by zero
end

%% Per-subcarrier MMSE equalization in TF domain
eqTF = zeros(N, M);
for col = 1:M
    H_col = hTF(:, col);       % Channel frequency response for this symbol
    % MMSE equalizer: W = H* / (|H|^2 + sigma^2)
    eqTF(:, col) = conj(H_col) ./ (abs(H_col).^2 + noiseVar) .* rxTF(:, col);
end

%% Convert equalized TF grid back to DD domain via SFFT
eqGrid = sqrt(M/N) * ifft(fft(eqTF, [], 1), [], 2);   % SFFT

%% Extract data symbols from equalized grid
dataSymbols = eqGrid(dataIdx);

end
