function [H, chInfo] = multipathChannel(cpSize, delta_f, inSig, velocity, ntnConfig)

%--------------------------------------------------------------------------
%
%   3GPP TR 38.811 NTN Channel Model with Shadowed-Rician LoS Fading
%
%   Uses MATLAB 5G Toolbox's nrTDLChannel in Custom mode with NTN-TDL
%   tap parameters from 3GPP TR 38.811 Section 6.9.2.
%
%   Compatible with all MATLAB versions that have nrTDLChannel.
%
%   Enhancement: LoS component uses Shadowed-Rician (Loo model) fading
%   to capture ionospheric scintillation and tropospheric shadowing.
%   Shadow std follows ITU-R P.681 elevation/scenario-dependent tables.
%
%   Supported profiles:
%     NTN-TDL-A  : NLOS, 3 taps  (Table 6.9.2-1)
%     NTN-TDL-C  : LOS,  derived from A via K-factor scaling (TR 38.901)
%     NTN-TDL-D  : LOS,  4 taps  (Table 6.9.2-4)
%
%--------------------------------------------------------------------------

%% Parse NTN configuration
if nargin < 5, ntnConfig = struct(); end
if ~isfield(ntnConfig, 'profile'),     ntnConfig.profile = 'NTN-TDL-C'; end
if ~isfield(ntnConfig, 'elevAngle'),   ntnConfig.elevAngle = 50; end
if ~isfield(ntnConfig, 'scenario'),    ntnConfig.scenario = 'Suburban'; end
if ~isfield(ntnConfig, 'fc'),          ntnConfig.fc = 2e9; end
if ~isfield(ntnConfig, 'altitude_km'), ntnConfig.altitude_km = 600; end

profile   = upper(ntnConfig.profile);
elevAngle = ntnConfig.elevAngle;
scenario  = ntnConfig.scenario;
fc        = ntnConfig.fc;

[N_sc, M_sym] = size(inSig);

%% ========================================================================
%  Physical parameters
%  ========================================================================
c_light    = physconst('LightSpeed');
v_ms       = velocity * 1e3 / 3600;
fd         = round(v_ms * fc / c_light);
sampleRate = N_sc * delta_f;
refDelaySpread_s = 30e-9;

%% ========================================================================
%  NTN-TDL Tap Tables (3GPP TR 38.811 Section 6.9.2)
%  ========================================================================

% NTN-TDL-A: NLOS, 3 taps (Table 6.9.2-1)
tdlA.delays_norm = [0,     1.0811, 2.8416];
tdlA.powers_dB   = [0,    -4.675, -6.482];

% NTN-TDL-D: LOS, 4 taps (Table 6.9.2-4)
tdlD.delays_norm = [0,     0.4740, 1.6470, 3.5090];
tdlD.powers_dB   = [0,    -7.127, -10.028, -13.590];

%% ========================================================================
%  Select profile, compute K-factor and tap parameters
%  ========================================================================
KFactor_dB = getNtnKFactor(elevAngle, scenario);

switch profile
    case 'NTN-TDL-A'
        baseDelays_s = tdlA.delays_norm * refDelaySpread_s;
        basePowers_dB = tdlA.powers_dB;
        isLOS = false;
        KFactor_dB = -Inf;

    case 'NTN-TDL-C'
        % LOS profile derived from NTN-TDL-A via K-factor scaling
        % (TR 38.901 Section 7.7.6)
        baseDelays_s = tdlA.delays_norm * refDelaySpread_s;
        basePowers_dB = tdlA.powers_dB;
        isLOS = true;

    case 'NTN-TDL-D'
        baseDelays_s = tdlD.delays_norm * refDelaySpread_s;
        basePowers_dB = tdlD.powers_dB;
        isLOS = true;

    otherwise
        error('Unsupported profile: %s. Use NTN-TDL-A, NTN-TDL-C, or NTN-TDL-D.', profile);
end

numPaths = length(baseDelays_s);

%% ========================================================================
%  K-factor scaling for LOS profiles
%  Scale all NLOS taps by 1/(K+1), LOS component K/(K+1) added separately
%  ========================================================================
if isLOS
    K_lin = 10^(KFactor_dB / 10);
    basePowers_lin = 10.^(basePowers_dB / 10);
    scatterPowers_lin = basePowers_lin / (K_lin + 1);
    scatterPowers_dB = 10 * log10(scatterPowers_lin);
else
    K_lin = 0;
    scatterPowers_dB = basePowers_dB;
end

%% ========================================================================
%  Generate fading via MATLAB nrTDLChannel (Custom mode, Rayleigh)
%  ========================================================================
channel = nrTDLChannel;
channel.DelayProfile         = 'Custom';
channel.PathDelays            = baseDelays_s;
channel.AveragePathGains      = scatterPowers_dB;
channel.FadingDistribution    = 'Rayleigh';
channel.MaximumDopplerShift   = 1;
channel.SampleRate            = sampleRate;
channel.NumTransmitAntennas   = 1;
channel.NumReceiveAntennas    = 1;

% Generate fading realization
numSamp = N_sc;
dummyIn = complex(randn(numSamp, 1), randn(numSamp, 1)) / sqrt(2);
[~, pathGainsAll] = channel(dummyIn);
% pathGainsAll: [numSamp x numPaths x 1 x 1]

% Extract base scatter gains at first sample
pathGains = squeeze(pathGainsAll(1, :, 1, 1)).';   % numPaths x 1

% Add Shadowed-Rician LOS component to first tap for LOS profiles
if isLOS
    % Shadowed-Rician (Loo model): LoS amplitude follows lognormal fading
    % Models ionospheric scintillation and tropospheric shadowing on the
    % direct path. Shadow std depends on elevation and environment.
    sigma_s_dB = getShadowingStd(elevAngle, scenario);

    % Lognormal shadow: 10^(N(0, sigma_s_dB) / 20)
    shadowFading_dB = sigma_s_dB * randn();
    shadowGain = 10^(shadowFading_dB / 20);

    % LoS amplitude with shadow fading + random phase
    losAmplitude = sqrt(K_lin / (K_lin + 1)) * shadowGain;
    losPhase = 2 * pi * rand();
    pathGains(1) = losAmplitude * exp(1j * losPhase) + pathGains(1);
end

%% ========================================================================
%  LEO Satellite Doppler Model
%  ========================================================================
scatterSpread = 0.05 * fd;

Vi = zeros(1, numPaths);
for i = 0:numPaths-1
    if numPaths > 1
        Vi(i+1) = fd + scatterSpread * cos(2 * pi * i / (numPaths - 1));
    else
        Vi(i+1) = fd;
    end
end

%% ========================================================================
%  Build TF-domain channel matrix H(subcarrier, symbol)
%  ========================================================================
T  = 1 / delta_f;
Ts = (1 + cpSize) / delta_f;

sc_idx  = (1:N_sc)' + N_sc/2;
sym_idx = 1:M_sym;

H = zeros(N_sc, M_sym);
for x = 1:numPaths
    hiPrime = pathGains(x) * (1 + 1i * pi * Vi(x) * T);
    expTerm = -2i * pi * (sc_idx * (delta_f * baseDelays_s(x)) ...
              - Vi(x) * Ts * sym_idx);
    H = H + hiPrime * exp(expTerm);
end

%% ========================================================================
%  Return channel parameters for DD-domain processing
%  ========================================================================
chInfo.pathDelays_s     = baseDelays_s;
chInfo.pathDopplers_Hz  = Vi;
chInfo.pathGains        = pathGains;
chInfo.numPaths         = numPaths;
chInfo.fd               = fd;
chInfo.KFactor_dB       = KFactor_dB;
chInfo.profile          = profile;
chInfo.elevAngle        = elevAngle;
chInfo.scenario         = scenario;
chInfo.refDelaySpread_s = refDelaySpread_s;

end


%% ========================================================================
%  LOCAL FUNCTION: Elevation-Angle-Dependent K-Factor
%  3GPP TR 38.811 Section 6.7.2, S-band LOS, per scenario
%  ========================================================================
function K_dB = getNtnKFactor(elevAngle, scenario)

elev_table = [10, 20, 30, 40, 50, 60, 70, 80, 90];

K_DenseUrban_dB = [ 4.4,  5.0,  5.4,  5.8,  6.2,  6.4,  6.6,  6.7,  6.8];
K_Urban_dB      = [ 7.0,  7.8,  8.3,  8.8,  9.2,  9.5,  9.7,  9.9, 10.0];
K_Suburban_dB   = [ 9.0, 10.0, 10.8, 11.4, 12.0, 12.4, 12.7, 12.9, 13.0];
K_Rural_dB      = [12.0, 13.2, 14.0, 14.6, 15.0, 15.4, 15.6, 15.8, 16.0];

switch scenario
    case 'DenseUrban'
        K_table = K_DenseUrban_dB;
    case 'Urban'
        K_table = K_Urban_dB;
    case 'Suburban'
        K_table = K_Suburban_dB;
    case 'Rural'
        K_table = K_Rural_dB;
    otherwise
        warning('Unknown scenario "%s", using Urban.', scenario);
        K_table = K_Urban_dB;
end

elevAngle = max(10, min(90, elevAngle));
K_dB = interp1(elev_table, K_table, elevAngle, 'pchip');

end


%% ========================================================================
%  LOCAL FUNCTION: Elevation-Angle-Dependent Shadow Fading Std
%  ITU-R P.681 / Loo model parameters for LEO satellite channels
%  ========================================================================
function sigma_dB = getShadowingStd(elevAngle, scenario)
% Shadow fading standard deviation for LoS path (dB)
% Based on ITU-R P.681 and Loo model parameters for LEO satellite channels.
% Higher elevation -> less shadowing; Rural -> less shadowing than Urban.

elev_table = [10, 20, 30, 40, 50, 60, 70, 80, 90];

sigma_DenseUrban = [4.0, 3.5, 3.0, 2.7, 2.4, 2.2, 2.0, 1.9, 1.8];
sigma_Urban      = [3.2, 2.8, 2.4, 2.1, 1.8, 1.6, 1.4, 1.3, 1.2];
sigma_Suburban    = [2.5, 2.1, 1.8, 1.5, 1.3, 1.1, 0.9, 0.8, 0.7];
sigma_Rural       = [1.5, 1.2, 1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.3];

switch scenario
    case 'DenseUrban'
        sigma_table = sigma_DenseUrban;
    case 'Urban'
        sigma_table = sigma_Urban;
    case 'Suburban'
        sigma_table = sigma_Suburban;
    case 'Rural'
        sigma_table = sigma_Rural;
    otherwise
        sigma_table = sigma_Urban;
end

elevAngle = max(10, min(90, elevAngle));
sigma_dB = interp1(elev_table, sigma_table, elevAngle, 'pchip');
end
