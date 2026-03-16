function [eqDataSyms, marginalPMF] = otfsEqualizerMP(y, H, dataIdx, alphabet, noiseVar, mpConfig)

%--------------------------------------------------------------------------
%
%   Memory-efficient Message Passing (MP) equalizer for OTFS (DD domain).
%
%   Implements Algorithm 1 from:
%     Raviteja et al., "Interference Cancellation and Iterative Detection
%     for Orthogonal Time Frequency Space (OTFS) Modulation," IEEE TWC 2018.
%
%   Memory optimization: uses per-variable-node belief (numData x Q)
%   instead of per-edge PMF (numEdges x Q), reducing peak memory from
%   ~500 MB to ~50 MB for typical OTFS grid sizes.
%
%   The extrinsic property is maintained in Step A (interference mean/var
%   subtraction) via the "compute total, subtract self" trick.
%   Step B accumulates likelihoods per-symbol (loop over Q) to avoid
%   creating any numEdges x Q intermediate matrix.
%
%--------------------------------------------------------------------------
% Input arguments:
%
% y          NM x 1 received signal vector (after pilot subtraction,
%            column-major vectorization of N x M DD grid)
% H          NM x NM sparse DD channel matrix (from buildDDChannelMatrix)
% dataIdx    numData x 1 linear indices of data positions in NM grid
% alphabet   Q x 1 complex QAM constellation points
% noiseVar   Noise variance (scalar, total power per DD element)
% mpConfig   (optional) Struct with fields:
%              .maxIter  - Maximum MP iterations (default 10)
%              .damping  - Damping factor 0 < delta < 1 (default 0.7)
%
%--------------------------------------------------------------------------
% Function returns:
%
% eqDataSyms   numData x 1 estimated data symbols (soft decision)
% marginalPMF  numData x Q marginal probability over constellation
%
%--------------------------------------------------------------------------

%% Configuration
if nargin < 6, mpConfig = struct(); end
if ~isfield(mpConfig, 'maxIter'), mpConfig.maxIter = 10; end
if ~isfield(mpConfig, 'damping'), mpConfig.damping = 0.7; end

maxIter = mpConfig.maxIter;
damping = mpConfig.damping;

NM = size(H, 1);
alphabet = alphabet(:);     % Q x 1 column vector
Q = length(alphabet);
numData = length(dataIdx);

%% Reduce system to data columns only
% After pilot subtraction, x is non-zero only at dataIdx
H_data = H(:, dataIdx);     % NM x numData sparse

%% Extract edge information from sparse matrix
[d_edges, c_edges, h_edges] = find(H_data);
% d_edges:  observation node indices (1..NM), numEdges x 1
% c_edges:  variable node indices (1..numData), numEdges x 1
% h_edges:  channel coefficients, numEdges x 1
numEdges = length(d_edges);

if numEdges == 0
    eqDataSyms = zeros(numData, 1);
    marginalPMF = ones(numData, Q) / Q;
    return;
end

%% Pre-compute
alpha_pow = abs(alphabet).^2;          % Q x 1
noiseVar = max(noiseVar, 1e-12);
y_d = y(d_edges);                      % numEdges x 1 (observations at edges)
h_abs2 = abs(h_edges).^2;             % numEdges x 1

%% Initialize per-node belief: uniform over alphabet
belief = ones(numData, Q) / Q;         % numData x Q  (~2 MB)

%% ===================== Iterative Message Passing =====================
for iter = 1:maxIter

    %% ------ Step A: Compute extrinsic mean & variance per edge ------
    % E[x_c] and Var[x_c] from current per-node belief
    E_x_node = belief * alphabet;              % numData x 1
    E_x2_node = belief * alpha_pow;            % numData x 1
    Var_x_node = E_x2_node - abs(E_x_node).^2;% numData x 1
    Var_x_node = max(Var_x_node, 0);

    % Map node values to edges
    E_x_edge = E_x_node(c_edges);              % numEdges x 1
    Var_x_edge = Var_x_node(c_edges);          % numEdges x 1

    % Weighted contributions per edge
    w_mean = h_edges .* E_x_edge;              % numEdges x 1
    w_var  = h_abs2 .* Var_x_edge;             % numEdges x 1

    % Total per observation node
    total_mean = accumarray(d_edges, w_mean, [NM, 1], @sum, 0);
    total_var  = accumarray(d_edges, w_var,  [NM, 1], @sum, 0) + noiseVar;

    % Extrinsic: subtract self-contribution (preserves BP extrinsic property)
    mu_dc     = total_mean(d_edges) - w_mean;          % numEdges x 1
    sigma2_dc = total_var(d_edges)  - w_var;            % numEdges x 1
    sigma2_dc = max(sigma2_dc, 1e-12);

    %% ------ Step B: Accumulate log-likelihoods per variable node ------
    % Process one constellation point at a time to avoid numEdges x Q matrix.
    % For each symbol a_j:
    %   residual_j = y[d] - mu_dc - h * a_j       (numEdges x 1)
    %   log_lik_j  = -|residual_j|^2 / sigma2_dc  (numEdges x 1)
    %   Accumulate into new_belief_log(c, j) via accumarray

    y_minus_mu = y_d - mu_dc;                          % numEdges x 1 (pre-compute once)
    new_belief_log = zeros(numData, Q);

    for j = 1:Q
        residual_j = y_minus_mu - h_edges * alphabet(j);       % numEdges x 1
        log_lik_j = -abs(residual_j).^2 ./ sigma2_dc;          % numEdges x 1
        new_belief_log(:, j) = accumarray(c_edges, log_lik_j, [numData, 1], @sum, 0);
    end

    % Log-Sum-Exp normalization for numerical stability
    log_max = max(new_belief_log, [], 2);                       % numData x 1
    new_belief = exp(new_belief_log - log_max);                 % numData x Q
    bsum = sum(new_belief, 2);
    bsum(bsum == 0) = 1;
    new_belief = new_belief ./ bsum;

    % Damping update (eq. 32, Raviteja 2018)
    belief = damping * new_belief + (1 - damping) * belief;

    % Re-normalize after damping
    belief = belief ./ sum(belief, 2);
end

%% ===================== Output =====================
marginalPMF = belief;

% Soft decision: expected symbol value from marginal
eqDataSyms = marginalPMF * alphabet;    % numData x 1

end
