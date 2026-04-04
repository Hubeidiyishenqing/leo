function eqTF = applyOFDM_ICICancel(rxTF, chInfo, scs, cpSize, fd_precomp, N, M, noiseVar, numIter)
% OFDM receiver with successive ICI cancellation (SIC).
%
% Uses time-domain reconstruction for ICI estimation:
%   Iter 0: Standard one-tap MMSE equalization
%   Iter i:
%     1. Take equalized symbols from previous iteration as TX estimate
%     2. Pass through full time-domain channel (applyChannelTF_ICI)
%        to get the full received signal estimate (including ICI)
%     3. Subtract the desired diagonal component: ICI = Y_hat - H_diag*X_hat
%     4. Clean received signal: Y_clean = Y_rx - ICI
%     5. Re-equalize with one-tap MMSE
%
% This approach is provably correct because applyChannelTF_ICI computes
% the exact time-domain per-sample Doppler phase rotation and FFT.

if nargin < 9, numIter = 2; end

T  = 1 / scs;
Ts = (1 + cpSize) / scs;
delta_tau = 1 / (N * scs);

n_idx = (0:N-1)';
m_idx = 0:M-1;

% Pre-compute per-path parameters
numP = chInfo.numPaths;
nu_res = zeros(numP, 1);
li_arr = zeros(numP, 1);
for i = 1:numP
    li_arr(i) = round(chInfo.pathDelays_s(i) / delta_tau);
    nu_res(i) = chInfo.pathDopplers_Hz(i) - fd_precomp;
end

% One-tap channel (diagonal) for MMSE — same as buildTF_OFDMeq
H_diag = zeros(N, M);
for i = 1:numP
    sincFactor = sinc(nu_res(i) * T) * exp(1j * pi * nu_res(i) * T);
    delayPhase = exp(-1j * 2 * pi * n_idx * li_arr(i) / N);
    dopplerPhase = exp(1j * 2 * pi * nu_res(i) * m_idx * Ts);
    H_diag = H_diag + chInfo.pathGains(i) * sincFactor * (delayPhase * dopplerPhase);
end

% Initial MMSE equalization (iteration 0)
eqTF = conj(H_diag) ./ (abs(H_diag).^2 + noiseVar) .* rxTF;

% Iterative SIC via time-domain reconstruction
for iter = 1:numIter
    % Step 1: Use previous equalized output as TX symbol estimate
    x_hat = eqTF;

    % Step 2: Pass x_hat through the full time-domain channel
    % This produces Y_hat = H_full * x_hat (including all ICI)
    Y_hat = applyChannelTF_ICI(x_hat, chInfo, scs, cpSize, fd_precomp);

    % Step 3: Estimate ICI = full channel output - diagonal-only component
    % Diagonal component: H_diag .* x_hat (the desired signal part)
    ici_est = Y_hat - H_diag .* x_hat;

    % Step 4: Subtract ICI from actual received signal
    cleaned = rxTF - ici_est;

    % Step 5: Re-equalize the cleaned signal
    eqTF = conj(H_diag) ./ (abs(H_diag).^2 + noiseVar) .* cleaned;
end

end
