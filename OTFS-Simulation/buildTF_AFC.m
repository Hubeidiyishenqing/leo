function H_afc = buildTF_AFC(chInfo, scs, cpSize, N, M, fd_precomp)
% Build TF channel matrix after AFC (Automatic Frequency Control)
% Used for MMSE equalization. The AFC removes bulk Doppler fd_precomp,
% leaving only residual scatter Doppler which causes ICI.
T  = 1 / scs;
Ts = (1 + cpSize) / scs;
sc_idx  = (1:N)' + N/2;
sym_idx = 1:M;
H_afc = zeros(N, M);
for x = 1:chInfo.numPaths
    Vi_res = chInfo.pathDopplers_Hz(x) - fd_precomp;
    hi = chInfo.pathGains(x) * (1 + 1j*pi*Vi_res*T);
    exT = -2j*pi*(sc_idx*(scs*chInfo.pathDelays_s(x)) - Vi_res*Ts*sym_idx);
    H_afc = H_afc + hi * exp(exT);
end
end
