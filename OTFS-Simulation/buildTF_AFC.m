function H_eq = buildTF_AFC(chInfo, scs, cpSize, N, M, fd_precomp)
% Build TF equalization channel matched to applyChannelDD integer mode.
%
% applyChannelDD (integer mode) applies 2D circular shifts in the DD
% domain: rxDD = sum_i h_i * circshift(txDD, [li, ki]).  For OFDM,
% the effective TF channel G = ISFFT * H_DD * SFFT is diagonal when
% H_DD is a 2D circulant.  The diagonal elements are:
%
%   G[n,m] = sum_i h_i * exp(j2pi*(n-1)*li/N) * exp(-j2pi*(m-1)*ki/M)
%
% where li = round(delay / delta_tau), ki = round((doppler-fd) / delta_nu).
%
% This function computes exactly these diagonal elements so that the
% TF-MMSE equalizer is perfectly consistent with the DD channel.

Ts = (1 + cpSize) / scs;          % OFDM symbol duration with CP
delta_tau = 1 / (N * scs);        % Delay resolution (s/bin)
delta_nu  = 1 / (M * Ts);         % Doppler resolution (Hz/bin)

n_idx = (0:N-1)';                 % 0-indexed subcarrier
m_idx = 0:M-1;                    % 0-indexed OFDM symbol

H_eq = zeros(N, M);
for i = 1:chInfo.numPaths
    li = round(chInfo.pathDelays_s(i) / delta_tau);
    ki = round((chInfo.pathDopplers_Hz(i) - fd_precomp) / delta_nu);
    hi = chInfo.pathGains(i);
    H_eq = H_eq + hi * (exp(1j*2*pi*n_idx*li/N) * exp(-1j*2*pi*m_idx*ki/M));
end
end
