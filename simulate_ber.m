function [ber,PMR] = simulate_ber(K, n, L, snr, num_tests)
    ber = zeros(1, K*n);
    N0 = (10^(-snr/10));
    sigma = sqrt(N0/2);
    
    
    for bit = 1:K*n
        wrong_dec = 0;
        for tests = 1:num_tests

            % generate inf. w   ord
            iwd = randi([0 1], K, n);
            %предположим, заморозим первые 10 битов
            iwd(1:10) = 0;

            % encode
            cwd = polar_transform(iwd);

            % modulate ans calculate sums
            tx = sum(1-2*cwd, 1);

            % add noise
            noise_vector = randn(1,n)*sigma;
            rx = tx + noise_vector;

            % demodulate
            temp_sums = sum(1-2*de2bi(0:2^K-1, K), 2);
            temp_sums = repmat(temp_sums, 1, n);
            in_logProbs = -abs(repmat(rx, 2^K, 1) - temp_sums).^2/2/sigma^2;                        

            frozen_matrix = 2*ones(K, n);
            %скажем, что мы знаем, что первые 10 битов точно равны 0
            frozen_matrix(1:10) = 0;
            % здесь если элемент матрицы равен 0, то мы точно знаем что там
            % 0, если 1, то мы точно знаем что там 1, если 2, то мы не
            % знаем что там и эт надо декодить. По сути так в декодер
            % передаётся информация о расстановке замороденых битов
%              frozen_matrix(1:bit-1) = iwd(1:bit-1);
            
            % decode
            [est_iwd, prob] = pcscl_jsc(K, log2(n), L, frozen_matrix, in_logProbs);
            
            
            
            est_iwd = transpose(de2bi(est_iwd, K));
            PMR = max(abs(prob))-min(abs(prob));
%             plot(PMR);

            est_iwddd = reshape(est_iwd,[],n);
            imshow([est_iwddd;iwd*2],[0 2])


            if (est_iwd(bit) ~= iwd(bit))
                wrong_dec = wrong_dec + 1;
            end
        end
        ber(bit) = wrong_dec/tests;
        fprintf('\tber (%d) = %f\n', bit, ber(bit));
    end
%     save(sprintf('ber_K=%d_n=%d_snr=%1.1f.mat', K, n, snr)); 
end
