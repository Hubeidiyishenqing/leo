function rxDD = applyChannelDD(txDD, chInfo, scs, cpSize, dopplerPrecomp_Hz, useFracDoppler)

%--------------------------------------------------------------------------
%
%   Applies the multipath channel directly in the Delay-Doppler domain
%   via 2D circular convolution.
%
%   Supports two modes:
%     Integer mode (default): rounds delay & Doppler to nearest bin.
%       Exact match with buildDDChannelMatrix(Ni=0). Use for BER evaluation.
%     Fractional mode (useFracDoppler=true): models both Inter-Delay
%       Interference and Inter-Doppler Interference via 2D Dirichlet kernel
%       spreading. Enables sub-bin Quinn estimation for both delay and
%       Doppler. Use for navigation accuracy evaluation.
%
%   Supports Doppler pre-compensation: subtracts a bulk Doppler shift
%   (e.g. from satellite ephemeris) so the residual Doppler is small
%   enough for the DD-domain pilot guard band to capture.
%
%--------------------------------------------------------------------------
% Input arguments:
%
% txDD              N x M transmit DD grid
% chInfo            Struct from multipathChannel with fields:
%                     .pathDelays_s    - path delays in seconds
%                     .pathDopplers_Hz - path Doppler shifts in Hz
%                     .pathGains       - path complex gains
%                     .numPaths        - number of channel paths
% scs               Subcarrier spacing (Hz)
% cpSize            Cyclic prefix ratio
% dopplerPrecomp_Hz (optional) Bulk Doppler to pre-compensate (Hz)
%                   Typically from satellite ephemeris. Default = 0.
% useFracDoppler    (optional) Enable fractional Doppler IDI model.
%                   Default = false (integer Doppler only).
%
%--------------------------------------------------------------------------
% Function returns:
%
% rxDD              N x M received DD grid after channel (pre-compensated)
%
%--------------------------------------------------------------------------

if nargin < 5, dopplerPrecomp_Hz = 0; end
if nargin < 6, useFracDoppler = false; end

[N, M] = size(txDD);

% DD grid resolutions
Ts = (1 + cpSize) / scs;           % OFDM symbol duration with CP (s)
delta_tau = 1 / (N * scs);         % Delay resolution (s/bin)
delta_nu  = 1 / (M * Ts);          % Doppler resolution (Hz/bin)

rxDD = zeros(N, M);

if ~useFracDoppler
    %% Integer Doppler mode (matches buildDDChannelMatrix with Ni=0)
    for i = 1:chInfo.numPaths
        li = round(chInfo.pathDelays_s(i) / delta_tau);
        ki = round((chInfo.pathDopplers_Hz(i) - dopplerPrecomp_Hz) / delta_nu);
        hi = chInfo.pathGains(i);
        rxDD = rxDD + hi * circshift(txDD, [li, ki]);
    end
else
    %% Fractional Delay + Fractional Doppler mode (2D spreading)
    %  Both delay and Doppler use Dirichlet kernel spreading coefficients.
    %  Delay spreading enables sub-bin delay estimation via Quinn estimator.
    %  Doppler spreading models IDI for sub-bin Doppler estimation.
    Ni_k = 2;    % Doppler spreading half-width
    Ni_l = 2;    % Delay spreading half-width
    l_idx = (0:N-1)';

    for i = 1:chInfo.numPaths
        hi = chInfo.pathGains(i);

        % Fractional delay decomposition
        li_frac = chInfo.pathDelays_s(i) / delta_tau;
        li = floor(li_frac);
        lambda = li_frac - li;       % fractional part [0, 1)

        % Fractional Doppler decomposition
        ki_frac = (chInfo.pathDopplers_Hz(i) - dopplerPrecomp_Hz) / delta_nu;
        ki = floor(ki_frac);
        kappa = ki_frac - ki;        % fractional part [0, 1)

        % Phase correction for fractional Doppler
        phase_col = exp(1j * 2 * pi * l_idx * (ki + kappa) / (N * M));

        % Check if both are near-integer (skip spreading)
        nearIntDelay = (abs(lambda) < 1e-10) || (abs(lambda - 1) < 1e-10);
        nearIntDoppler = (abs(kappa) < 1e-10) || (abs(kappa - 1) < 1e-10);

        if nearIntDelay && nearIntDoppler
            li_int = round(li_frac);
            ki_int = round(ki_frac);
            rxDD = rxDD + hi * circshift(txDD, [li_int, ki_int]);
        else
            for p = -Ni_l:Ni_l
                % Delay spreading coefficient (Dirichlet kernel)
                if abs(lambda - p) < 1e-10
                    alpha_p = 1;
                else
                    alpha_p = (1/N) * (1 - exp(1j*2*pi*(lambda - p))) / ...
                              (1 - exp(1j*2*pi*(lambda - p)/N));
                end
                if abs(alpha_p) < 1e-15, continue; end

                for q = -Ni_k:Ni_k
                    % Doppler spreading coefficient (Dirichlet kernel)
                    if abs(kappa - q) < 1e-10
                        beta_q = 1;
                    else
                        beta_q = (1/M) * (1 - exp(1j*2*pi*(kappa - q))) / ...
                                 (1 - exp(1j*2*pi*(kappa - q)/M));
                    end
                    if abs(beta_q) < 1e-15, continue; end

                    shifted = circshift(txDD, [li + p, ki + q]);
                    rxDD = rxDD + hi * alpha_p * beta_q * (phase_col .* shifted);
                end
            end
        end
    end
end

end
