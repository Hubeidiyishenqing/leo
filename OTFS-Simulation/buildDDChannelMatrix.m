function H = buildDDChannelMatrix(pathDelays, pathDopplers, pathGains, N, M, Ni)

%--------------------------------------------------------------------------
%
%   Builds the sparse DD-domain channel matrix H (NM x NM)
%   for the input-output relation y = Hx + z.
%
%   Supports integer Doppler (Ni=0, matches applyChannelDD.m) and
%   fractional Doppler with IDI truncation (Ni>0, per Raviteja 2018).
%
%   The DD-domain channel acts via 2D circular convolution:
%     y[l',k'] = sum_i h_i * x[(l'-l_i) mod N, (k'-k_i) mod M]
%
%   For fractional Doppler (Ni>0), IDI spreading adds extra taps:
%     y[l',k'] = sum_i h_i * sum_{q=-Ni}^{Ni} beta_i(q) * phase_i *
%                x[(l'-l_i) mod N, (k'-k_i-q) mod M]
%
%--------------------------------------------------------------------------
% Input arguments:
%
% pathDelays     L x 1, delay shifts in bins (integer)
% pathDopplers   L x 1, Doppler shifts in bins (integer for Ni=0,
%                or full fractional value for Ni>0)
% pathGains      L x 1, complex path gains
% N              Number of delay bins (rows of DD grid)
% M              Number of Doppler bins (columns of DD grid)
% Ni             (optional) IDI truncation parameter, default = 0
%                Ni=0: integer Doppler only (fast, matches applyChannelDD)
%                Ni>0: fractional Doppler with 2*Ni+1 spread taps
%
%--------------------------------------------------------------------------
% Function returns:
%
% H              NM x NM sparse channel matrix
%
%--------------------------------------------------------------------------

if nargin < 6
    Ni = 0;
end

NM = N * M;
numPaths = length(pathGains);

% Create 0-based 2D index grids (vectorized over all NM positions)
[l_grid, k_grid] = ndgrid(0:N-1, 0:M-1);
l_vec = l_grid(:);   % NM x 1, delay indices (0-based)
k_vec = k_grid(:);   % NM x 1, Doppler indices (0-based)
c_vec = (1:NM)';     % Input (column) indices (1-based)

% Estimate max non-zeros and pre-allocate
numTaps = 2 * Ni + 1;
maxNNZ = numPaths * numTaps * NM;
rowAll = zeros(maxNNZ, 1);
colAll = zeros(maxNNZ, 1);
valAll = zeros(maxNNZ, 1);
count = 0;

for i = 1:numPaths
    hi = pathGains(i);

    if Ni == 0
        % --- Integer Doppler: simple circular shift ---
        li = round(pathDelays(i));
        ki = round(pathDopplers(i));

        l_out = mod(l_vec + li, N);
        k_out = mod(k_vec + ki, M);
        d_vec = l_out + k_out * N + 1;   % 1-based output indices

        idx = count + (1:NM);
        rowAll(idx) = d_vec;
        colAll(idx) = c_vec;
        valAll(idx) = hi;
        count = count + NM;
    else
        % --- Fractional Doppler: IDI spreading ---
        li = round(pathDelays(i));            % Integer delay (bins)
        ki_full = pathDopplers(i);            % Full Doppler (possibly fractional)
        ki = floor(ki_full);                  % Integer Doppler part
        kappa_i = ki_full - ki;               % Fractional Doppler part [0, 1)

        for q = -Ni:Ni
            % IDI coefficient (sinc-like spreading function)
            if abs(kappa_i - q) < 1e-10
                beta_q = 1;
            else
                beta_q = (1/M) * (1 - exp(1j*2*pi*(kappa_i - q))) / ...
                         (1 - exp(1j*2*pi*(kappa_i - q)/M));
            end

            coeff_base = hi * beta_q;
            if abs(coeff_base) < 1e-15
                continue;
            end

            % Output positions via circular shift
            l_out = mod(l_vec + li, N);
            k_out = mod(k_vec + ki + q, M);
            d_vec = l_out + k_out * N + 1;

            % Phase factor for rectangular pulse (eq. 24, Raviteja 2018)
            % phase = exp(j*2*pi*(l_out)*(ki + kappa_i) / (N*M))
            phase_vec = exp(1j * 2 * pi * l_out .* (ki + kappa_i) / (N * M));

            idx = count + (1:NM);
            rowAll(idx) = d_vec;
            colAll(idx) = c_vec;
            valAll(idx) = coeff_base * phase_vec;
            count = count + NM;
        end
    end
end

% Trim unused pre-allocation and build sparse matrix
rowAll = rowAll(1:count);
colAll = colAll(1:count);
valAll = valAll(1:count);
H = sparse(rowAll, colAll, valAll, NM, NM);

end
