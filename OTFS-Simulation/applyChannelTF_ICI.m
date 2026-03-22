function rxTF = applyChannelTF_ICI(txTF, chInfo, scs, cpSize, dopplerPrecomp_Hz)
% Apply multipath channel to OFDM signal with intra-symbol ICI modeling.
%
% Unlike applyChannelDD (which applies 2D circular shifts in DD domain,
% giving a perfectly diagonal TF channel), this function processes each
% OFDM symbol in the time domain with per-sample Doppler phase rotation.
% After FFT at the receiver, inter-carrier interference (ICI) from
% Doppler naturally appears as off-diagonal terms.
%
% Processing per OFDM symbol m:
%   1. IFFT: subcarrier data -> time-domain samples
%   2. Multipath + Doppler: each path has delay (circular shift) and
%      per-sample phase rotation exp(j*2*pi*nu_i*n*Ts_sample)
%   3. Receiver CFO correction: derotate by bulk Doppler estimate
%   4. FFT: time-domain -> received subcarriers (with residual ICI)
%
% The receiver CFO correction removes inter-symbol Doppler rotation
% but cannot undo intra-symbol ICI that occurred before the FFT.
%
% Inputs:
%   txTF              N x M transmit TF grid (subcarrier x OFDM symbol)
%   chInfo            Channel struct from multipathChannel
%   scs               Subcarrier spacing (Hz)
%   cpSize            CP ratio
%   dopplerPrecomp_Hz Bulk Doppler for receiver CFO correction (Hz)
%
% Output:
%   rxTF              N x M received TF grid (with ICI from Doppler)

if nargin < 5, dopplerPrecomp_Hz = 0; end

[N, M] = size(txTF);
Ts_sample = 1 / (N * scs);           % Sample period (s)
Ts_sym    = (1 + cpSize) / scs;      % OFDM symbol duration with CP (s)
delta_tau = 1 / (N * scs);           % Delay resolution (s/bin)

rxTF = zeros(N, M);
n_vec = (0:N-1)';                    % Sample indices within symbol

for m = 1:M
    % 1. IFFT to time domain
    x_td = ifft(txTF(:, m));

    % 2. Apply multipath channel with per-sample Doppler
    y_td = zeros(N, 1);
    for i = 1:chInfo.numPaths
        li = round(chInfo.pathDelays_s(i) / delta_tau);
        nu_i = chInfo.pathDopplers_Hz(i);  % Full Doppler (Hz)

        % Delay: circular shift in time domain
        x_delayed = circshift(x_td, li);

        % Doppler: per-sample phase rotation (causes ICI within symbol)
        % Plus inter-symbol phase at symbol start time
        t_samples = n_vec * Ts_sample + (m - 1) * Ts_sym;
        dopplerPhase = exp(1j * 2 * pi * nu_i * t_samples);

        y_td = y_td + chInfo.pathGains(i) * (dopplerPhase .* x_delayed);
    end

    % 3. Receiver CFO correction: remove bulk Doppler estimate
    %    This corrects inter-symbol rotation but NOT intra-symbol ICI
    t_corr = n_vec * Ts_sample + (m - 1) * Ts_sym;
    cfoCorrection = exp(-1j * 2 * pi * dopplerPrecomp_Hz * t_corr);
    y_corrected = y_td .* cfoCorrection;

    % 4. FFT back to frequency domain (ICI appears here)
    rxTF(:, m) = fft(y_corrected);
end

end
