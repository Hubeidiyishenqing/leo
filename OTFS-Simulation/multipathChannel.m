function [H, chInfo] = multipathChannel(cpSize, delta_f, inSig, velocity)

%--------------------------------------------------------------------------
%
%   Generates the time-frequency (TF) domain channel transfer function
%   for an LEO satellite Rician fading channel with Doppler shift.
%
%   Returns a 2D matrix H(subcarrier, symbol) for element-wise
%   multiplication with the TF-domain signal (after ISFFT, before OFDM mod).
%
%--------------------------------------------------------------------------
% Input arguments:
%
% cpSize                        OFDM cyclic prefix ratio
% delta_f                       OFDM subcarrier spacing (Hz)
% inSig                         Matrix whose size defines the TF grid (N_sc x M_sym)
% velocity                      Channel mobility in km/hr
%
%--------------------------------------------------------------------------
% Function returns:
%
% H                             N_sc x M_sym TF-domain channel matrix
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

% Get TF grid dimensions
[N_sc, M_sym] = size(inSig);       % subcarriers x OFDM symbols

% Generate Channel Parameters (LEO satellite Rician channel)
maxDelayspread = min(0.5*((cpSize)/delta_f), 500e-9); % Cap at 500 ns for LEO
L = min(round(2*maxDelayspread * N_sc*delta_f), 6);   % Cap at 6 paths for LEO
L = max(L, 2);                                         % At least 2 paths
step = maxDelayspread/L;                  % calculate difference between delays
pathDelays = (0:step:maxDelayspread);     % Discrete even delays of L-path channel
range shuffle;                            % Shuffle random no. generator

% Rician K-factor for LEO satellite (strong LoS component)
K_ric_dB = 10;                              % K-factor in dB
K_ric = 10^(K_ric_dB/10);                  % Linear K-factor
avgPathGains = zeros(L, 1);
avgPathGains(1) = sqrt(K_ric / (K_ric + 1)); % LoS path (dominant)
% NLOS paths: exponentially decaying power profile
nlosPowerTotal = 1 / (K_ric + 1);
nlosDecay = exp(-(0:L-2) / max(1, (L-2)/2));
nlosDecay = nlosDecay / sum(nlosDecay);       % Normalize to sum=1
for i = 2:L
    avgPathGains(i) = sqrt(nlosPowerTotal * nlosDecay(i-1));
end

% Calculate Max Doppler Shift
v = velocity*1e3/3600;                   % Mobile speed (m/s)
fc = 2e9;                                % Carrier frequency (S-band for LEO)
fd = round(v*fc/physconst('lightspeed'));% Maximum Doppler shift to nearest Hz

% Generate Doppler shifts for LEO satellite channel
% All paths share the bulk Doppler shift fd from satellite motion.
% Multipath scattering adds a small random perturbation on top.
scatterSpread = 0.05 * fd;   % scattering spread << fd (5% of fd)
Vi = zeros(1, L);
for l=0:L-1
    Vi(l+1) = fd + scatterSpread * cos( (2*pi*l)/(L-1) );
end

% Initialize channel variables
T = 1/delta_f;                  % unextended OFDM symbol period
Ts = (1+cpSize)/delta_f;        % OFDM symbol period with CP
Ti = pathDelays;                % Path delays
hi = avgPathGains;              % Path gains

% Return channel parameters for DD-domain processing
chInfo.pathDelays_s = Ti;       % Path delays in seconds
chInfo.pathDopplers_Hz = Vi;    % Path Doppler shifts in Hz
chInfo.pathGains = hi;          % Path complex gains
chInfo.numPaths = L;            % Number of paths
chInfo.fd = fd;                 % Max Doppler shift (Hz)

% Create TF-domain channel matrix H(subcarrier, symbol)
% H(k, n) = sum_i hi(i) * (1 + j*pi*Vi(i)*T) * exp(-j2pi*(k_freq*df*Ti(i) - Vi(i)*n*Ts))
% where k_freq is the centered subcarrier index
H = zeros(N_sc, M_sym);

for sc = 1:N_sc                 % Loop over subcarriers (rows)
    for sym = 1:M_sym           % Loop over OFDM symbols (columns)
        for x = 1:L             % Sum over L multipath components
            % Centered subcarrier frequency index
            freqIdx = sc + N_sc/2;
            % TF channel: delay causes frequency-dependent phase,
            %             Doppler causes time-dependent phase
            expTerm = (-2*1i*(pi)) * (freqIdx.*delta_f.*Ti(x) - Vi(x).*(sym).*Ts);
            hiPrime = hi(x)*(1 + 1i*(pi).*Vi(x).*T);
            H(sc, sym) = H(sc, sym) + exp(expTerm) * hiPrime;
        end
    end
end

end
