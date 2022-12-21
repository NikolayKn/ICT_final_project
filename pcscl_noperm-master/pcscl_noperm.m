function [decoded_bits, llr_out, path_metrics] = pcscl_noperm(llr, frozen_bit_pattern, list_size)
% Run a polar list decoder and return the whole candidate list of size L
% (decoded bits, LLR and path metrics)

% Path metric type, see https://arxiv.org/pdf/1411.7282.pdf for more
% details. If true, the approximate path metric type is used, see eq. (14)
% and (15) for more details. If false -- exact path metrics are used, see
% eq. (12) and (13).
path_metric_approx = true;
min_sum = true; % Use min-sum check node operations by default

if length(size(frozen_bit_pattern)) ~= 2 || size(frozen_bit_pattern, 1) ~= 1
    error('frozen_bit_pattern must be non-empty one-dimensional array');
end
n_polar = size(frozen_bit_pattern, 2);
if 2^round(log2(n_polar)) ~= n_polar
    error('The message length must be a power of 2 for polar codes.');
end

[decoded_bits, llr_out, path_metrics] = pcscl_impl(...
    log2(n_polar), list_size, path_metric_approx, min_sum, ...
    double(frozen_bit_pattern), llr...
    );
% Maximum list size is 2^32 or the codebook size
k = max(size(frozen_bit_pattern(~frozen_bit_pattern))); % Information bits
max_lists = 2^min(31, k);

if 	(max_lists < (list_size))
    decoded_bits = decoded_bits(1:max_lists,:);
    llr_out      = llr_out(1:max_lists,:);
    path_metrics = path_metrics(1:max_lists,:);
end

end

