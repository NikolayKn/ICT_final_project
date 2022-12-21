K = 1; % number of users
n = 128; % length of a codeword
L = 16; % number of lists
snr = 100; % signal to noise ratio
num_tests = 100; % number of iterations

ser = simulate_ber(K, n, L, snr, num_tests);
plot (ser)