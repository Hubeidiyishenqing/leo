function rxDD = applyChannelDD(txDD, chInfo, scs, cpSize, dopplerPrecomp_Hz, useFracDoppler)

%--------------------------------------------------------------------------
%
%   Applies the multipath channel directly in the Delay-Doppler domain
%   via 2D circular convolution.
%
%   Supports two modes:
%     Integer Doppler (default): rounds to nearest bin, exact match with
%       buildDDChannelMatrix(Ni=0). Use for BER evaluation.
%     Fractional Doppler (useFracDoppler=true): models IDI spreading,
%       enables accurate fractional Doppler estimation by Quinn estimator.
%       Use for navigation accuracy evaluation.
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
    %% Fractional Doppler mode (IDI spreading, Ni=2)
    Ni = 2;
    l_idx = (0:N-1)';

    for i = 1:chInfo.numPaths
        hi = chInfo.pathGains(i);
        li = round(chInfo.pathDelays_s(i) / delta_tau);
        ki_frac = (chInfo.pathDopplers_Hz(i) - dopplerPrecomp_Hz) / delta_nu;
        ki = floor(ki_frac);
        kappa = ki_frac - ki;

        if abs(kappa) < 1e-10
            rxDD = rxDD + hi * circshift(txDD, [li, ki]);
        else
            phase_col = exp(1j * 2 * pi * l_idx * (ki + kappa) / (N * M));
            for q = -Ni:Ni
                if abs(kappa - q) < 1e-10
                    beta_q = 1;
                else
                    beta_q = (1/M) * (1 - exp(1j*2*pi*(kappa - q))) / ...
                             (1 - exp(1j*2*pi*(kappa - q)/M));
                end
                if abs(beta_q) < 1e-15, continue; end
                shifted = circshift(txDD, [li, ki + q]);
                rxDD = rxDD + hi * beta_q * (phase_col .* shifted);
            end
        end
    end
end

end
