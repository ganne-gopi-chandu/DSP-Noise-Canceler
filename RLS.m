function snr_value = compute_snr(original_signal, noisy_signal)
    noise = noisy_signal - original_signal;
    signal_power = mean(original_signal.^2);
    noise_power = mean(noise.^2);
    snr_value = 10 * log10(signal_power / noise_power);
    fprintf('SNR: %.2f dB\n', snr_value);
end

s = load("clean_speech.txt"); 
w = load("external_noise.txt");
d = load("noisy_speech.txt"); 
Fs = 44100;

N = length(d);
L = 15;
lambda = 0.9999;
delta = 1e4; 
h = zeros(L, 1); 
P = delta * eye(L);

estimated_noise = zeros(N, 1); 
output_signal = zeros(N, 1);   
buffer = zeros(L, 1);

for n = 1:N
    buffer = [w(n); buffer(1:end-1)];
    
    % Compute gain vector K
    K = (P * buffer) / (lambda + buffer' * P * buffer);
    
    % Compute error signal e_k
    e_k = d(n) - h' * buffer;
    
    % Update filter coefficients h
    h = h + K * e_k;
    
    % Update error covariance matrix P
    P = (P - K * buffer' * P) / lambda;
    
    estimated_noise(n) = h' * buffer;
    output_signal(n) = d(n) - estimated_noise(n);
end

snr_initial = compute_snr(s, d);
snr_final = compute_snr(s, output_signal);

%sound(output_signal, Fs);
