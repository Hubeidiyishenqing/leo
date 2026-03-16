function [ddGrid, dataIdx, pilotIdx, guardIdx, pilotInfo] = pilotPatternDD(dataSymbols, N, M, pilotConfig)

%--------------------------------------------------------------------------
%
%   Places pilot, guard bands, and data symbols on the Delay-Doppler grid
%   for OTFS integrated communication and navigation (LEO satellite).
%
%   Uses the "Impulse Pilot + Guard Band" scheme in the DD domain.
%   The pilot is a single high-energy impulse at (kp, lp), surrounded
%   by a rectangular zero-value guard region. Data symbols fill the rest.
%
%--------------------------------------------------------------------------
% Input arguments:
%
% dataSymbols       Column vector of QAM data symbols to place on the grid
% N                 Number of rows (delay bins / subcarriers)
% M                 Number of columns (Doppler bins / OFDM symbols)
% pilotConfig       Struct with fields:
%   .kp             Pilot Doppler index (1-based), default = ceil(M/2)
%   .lp             Pilot delay index (1-based), default = ceil(N/2)
%   .kGuard         Doppler guard band width (one side), default = 3
%   .lGuard         Delay guard band width (one side), default = 2
%   .pilotBoostdB   Pilot power boost in dB relative to avg data power
%   .totalPowerdBm  (optional) Total transmit power constraint in dBm
%
%--------------------------------------------------------------------------
% Function returns:
%
% ddGrid            N x M matrix: the delay-Doppler grid with pilot+data
% dataIdx           Linear indices of data symbol positions
% pilotIdx          Linear index of the pilot position
% guardIdx          Linear indices of guard band (zero) positions
% pilotInfo         Struct with pilot/guard metadata for receiver
%
%--------------------------------------------------------------------------

%% Default pilot configuration
if nargin < 4 || isempty(pilotConfig)
    pilotConfig = struct();
end
if ~isfield(pilotConfig, 'kp'),           pilotConfig.kp = ceil(M/2); end
if ~isfield(pilotConfig, 'lp'),           pilotConfig.lp = ceil(N/2); end
if ~isfield(pilotConfig, 'kGuard'),       pilotConfig.kGuard = 3;     end
if ~isfield(pilotConfig, 'lGuard'),       pilotConfig.lGuard = 2;     end
if ~isfield(pilotConfig, 'pilotBoostdB'), pilotConfig.pilotBoostdB = 10; end

kp = pilotConfig.kp;
lp = pilotConfig.lp;
kGuard = pilotConfig.kGuard;
lGuard = pilotConfig.lGuard;

%% Validate guard band fits within the grid
if (kp - kGuard < 1) || (kp + kGuard > M)
    error('Doppler guard band exceeds grid boundary. Reduce kGuard or adjust kp.');
end
if (lp - lGuard < 1) || (lp + lGuard > N)
    error('Delay guard band exceeds grid boundary. Reduce lGuard or adjust lp.');
end

%% Initialize grid
ddGrid = zeros(N, M);

%% Identify guard region indices (rectangular block around pilot)
guardRows = (lp - lGuard):(lp + lGuard);   % delay indices
guardCols = (kp - kGuard):(kp + kGuard);   % Doppler indices

% Create logical mask for guard+pilot region
guardMask = false(N, M);
guardMask(guardRows, guardCols) = true;

% Pilot position
pilotMask = false(N, M);
pilotMask(lp, kp) = true;

% Guard = guard region minus the pilot itself
guardOnlyMask = guardMask & ~pilotMask;

% Data positions = everything outside the guard region
dataMask = ~guardMask;

%% Get linear indices
pilotIdx = find(pilotMask);
guardIdx = find(guardOnlyMask);
dataIdx  = find(dataMask);

%% Calculate number of available data positions
numDataPositions = length(dataIdx);

%% Place data symbols
% Truncate or zero-pad data to fit available positions
if length(dataSymbols) > numDataPositions
    % More data than available slots: truncate
    ddGrid(dataIdx) = dataSymbols(1:numDataPositions);
elseif length(dataSymbols) < numDataPositions
    % Fewer data symbols: place what we have, rest stays zero
    ddGrid(dataIdx(1:length(dataSymbols))) = dataSymbols;
else
    ddGrid(dataIdx) = dataSymbols;
end

%% Place pilot symbol with power boost
% Pilot amplitude = sqrt(boostLinear) * avg_data_amplitude
% This ensures pilot energy is boosted relative to data
avgDataPower = mean(abs(dataSymbols).^2);
pilotBoostLinear = 10^(pilotConfig.pilotBoostdB / 10);
pilotAmplitude = sqrt(pilotBoostLinear * avgDataPower);
ddGrid(lp, kp) = pilotAmplitude;   % Real-valued impulse pilot

%% Power normalization under total power constraint
% Total power = pilot power + data power
% Normalize so that total average power per symbol = 1
totalPower = sum(abs(ddGrid(:)).^2);
numActiveSymbols = numDataPositions + 1;  % data + pilot
if totalPower > 0
    normFactor = sqrt(numActiveSymbols / totalPower);
    ddGrid = ddGrid * normFactor;
end

%% Store pilot info for receiver side
pilotInfo.kp = kp;
pilotInfo.lp = lp;
pilotInfo.kGuard = kGuard;
pilotInfo.lGuard = lGuard;
pilotInfo.numDataPositions = numDataPositions;
pilotInfo.numGuardPositions = length(guardIdx);
pilotInfo.pilotBoostdB = pilotConfig.pilotBoostdB;
pilotInfo.overheadPercent = 100 * (length(guardIdx) + 1) / (N * M);
pilotInfo.pilotPower = abs(ddGrid(lp, kp))^2;
pilotInfo.avgDataPower = mean(abs(ddGrid(dataIdx)).^2);
pilotInfo.effectivePilotBoostdB = 10*log10(pilotInfo.pilotPower / pilotInfo.avgDataPower);

end
