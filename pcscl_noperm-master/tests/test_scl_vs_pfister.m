% Run successive cancellaion list and check that output LLR are exactly the
% same as the Pfister's decoder (succeiive cancellation)
% http://pfister.ee.duke.edu/courses/ecen655/polar.pdf
% The only change in Pfister's algorithm is the min-sum based check layer
% in LLR domain (see polar_decode_pfister function)

n_polar = 128;
frozen_pattern = randi([0, 1], 1, n_polar);
k = sum(frozen_pattern == 0);

iwd = randi([0, 1], 1, k);
iwd_polar = zeros(1, n_polar);
iwd_polar(frozen_pattern == 0) = iwd;
cwd = polar_transform_noperm(iwd_polar);
list_size = 1;
bps = 2;
sigma = 0.5;

% Run BPSK channel
tx = 1 - 2 * cwd;
rx = tx + randn(size(tx)) * sigma;
llr = 2 * rx / sigma^2;
% Use recursive SC decoder implementation
[llr_expected, cwd_expected] = polar_sc_recursive_llr_noperm(llr, frozen_pattern);
% Run SCL decoder
[decoded_bits, llr_got, ~] = pcscl_noperm(llr, frozen_pattern, list_size);
% Check codeword
assert(sum(polar_transform_noperm(cwd_expected) == decoded_bits) == n_polar);
% Check LLR for decoded bits
assert(max(abs(llr_got - llr_expected)) < 1e-3);
