% Run successive cancellaion list and check that output LLR are exactly the
% same as the Pfister's decoder (succeiive cancellation)
% http://pfister.ee.duke.edu/courses/ecen655/polar.pdf
% The only change in Pfister's algorithm is the min-sum based check layer
% in LLR domain (see polar_decode_pfister function)

n_polar = 32;
frozen_pattern = ones(1, n_polar);
k = 12; % The number of information bits
list_size = 2^k;
inf_idx = randsample(1:n_polar, k);
frozen_pattern(inf_idx) = 0;

iwd = randi([0, 1], 1, k);
iwd_polar = zeros(1, n_polar);
iwd_polar(frozen_pattern == 0) = iwd;
cwd = polar_transform_noperm(iwd_polar);

iwd_book = logical(de2bi(0:(2^k - 1), k));
iwd_book_polar = zeros(2^k, n_polar);
iwd_book_polar(:, frozen_pattern == 0) = iwd_book;
cwd_book = polar_transform_noperm(iwd_book_polar);

bps = 2;
sigma = 0.5;

% Run BPSK channel
tx = 1 - 2 * cwd;
rx = tx + randn(size(tx)) * sigma;
llr = 2 * rx / sigma^2;
% Use List decoder with the list size equal to codebook size
[decoded_bits, llr_got, path_metrics] = pcscl_noperm(llr, frozen_pattern, list_size);
[~, index] = max(path_metrics);
iwd_hat_polar = decoded_bits(index, :);
% Derive ML-based solution
scores = cwd_book * llr';
[~, id] = min(scores);
cwd_ml = cwd_book(id, :);
% Check that result is the same
assert(sum(cwd_ml == polar_transform_noperm(iwd_hat_polar)) == n_polar);