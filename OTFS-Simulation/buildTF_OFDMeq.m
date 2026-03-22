function H_eq = buildTF_OFDMeq(chInfo, scs, cpSize, N, M, fd_precomp)
% Build TF equalization channel for OFDM one-tap equalizer.
%
% This matches the main diagonal of the OFDM TF channel after receiver
% CFO correction. It includes the sinc amplitude loss from intra-symbol
% Doppler (ICI steals energy from the main tap) and phase rotation
% across subcarriers and symbols.
%
% Delay is rounded to integer samples to match applyChannelTF_ICI (which
% uses circshift for delay). Doppler remains continuous (per-sample phase
% rotation in applyChannelTF_ICI is inherently continuous).
%
% H_eq[n,m] = sum_i h_i * sinc(nu_res_i * T) * exp(j*pi*nu_res_i*T)
%             * exp(-j*2*pi*n*li/N) * exp(j*2*pi*nu_res_i*m*Ts)
%
% Inputs:
%   chInfo        Channel struct with pathDelays_s, pathDopplers_Hz, pathGains
%   scs           Subcarrier spacing (Hz)
%   cpSize        CP ratio
%   N             Number of subcarriers
%   M             Number of OFDM symbols
%   fd_precomp    Bulk Doppler for CFO correction (Hz)

T  = 1 / scs;                        % Useful OFDM symbol duration
Ts = (1 + cpSize) / scs;             % Full symbol duration with CP
delta_tau = 1 / (N * scs);           % Delay resolution (s/sample)

n_idx = (0:N-1)';                    % Subcarrier indices
m_idx = (0:M-1);                     % Symbol indices

H_eq = zeros(N, M);
for i = 1:chInfo.numPaths
    hi    = chInfo.pathGains(i);
    tau_i = chInfo.pathDelays_s(i);
    nu_res = chInfo.pathDopplers_Hz(i) - fd_precomp;  % Residual Doppler

    % Sinc loss from intra-symbol Doppler (ICI energy leakage)
    sincFactor = sinc(nu_res * T) * exp(1j * pi * nu_res * T);

    % Subcarrier phase from delay (integer samples, matching applyChannelTF_ICI)
    li = round(tau_i / delta_tau);
    delayPhase = exp(-1j * 2 * pi * n_idx * li / N);

    % Symbol phase from residual Doppler (continuous)
    dopplerPhase = exp(1j * 2 * pi * nu_res * m_idx * Ts);

    H_eq = H_eq + hi * sincFactor * (delayPhase * dopplerPhase);
end

end
