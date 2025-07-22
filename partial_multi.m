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
lambda = 1;
delta = 1e4; 
h = zeros(L, 1); 
P = delta * eye(L);

estimated_noise = zeros(N, 1); 
filtered_noise = zeros(N, 1);
output_signal_partial = zeros(N, 1);  
output_signal_full=zeros(N,1);
buffer = zeros(L, 1);

% Notch filter parameters
notch_freqs = [];%[1000,2734.3,1479.06];
r = 0.999; 

Nfreq = length(notch_freqs);

notch_x = zeros(3, 1);        
notch_y_all = zeros(3, Nfreq);

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
    
    % Get the estimated noise
    estimated_noise(n) = h' * buffer;
    
    % Apply notch filter to the estimated noise

    [tonal_sum,notch_x,notch_y_all] = multi_notch_filter_step(estimated_noise(n), notch_freqs, Fs, r, notch_x, notch_y_all);

    filtered_noise(n) = tonal_sum;

    % Subtract tonal components from estimated noise
    output_signal_partial(n) = d(n) - estimated_noise(n)+ filtered_noise(n);
    output_signal_full(n) = d(n) - estimated_noise(n);

end

function [tonal_sum,notch_x,notch_y_all] = multi_notch_filter_step(input_sample, notch_freqs, Fs, r, notch_x, notch_y_all)
    
    Nfreq = length(notch_freqs);
    tonal_components = zeros(Nfreq, 1);

    notch_x = [input_sample; notch_x(1:2)];

    for k = 1:Nfreq
        notch_y = notch_y_all(:, k);
        notch_y(2:3) = notch_y(1:2);

        a = 2 * cos(2 * pi * (notch_freqs(k) / Fs));
        b = [1, -a, 1];
        a_coeff = [1, -a*r, r^2];

        tonal_components(k) = b * notch_x - a_coeff(2:end) * notch_y(2:3);

        notch_y(1) = tonal_components(k);
        notch_y_all(:, k) = notch_y;
        tonal_components(k)= input_sample- tonal_components(k);
    end

    tonal_sum = sum(tonal_components);
end

snr_initial = compute_snr(s, d);
snr_final = compute_snr(s, output_signal_partial);

% Uncomment to listen to the result
 sound(output_signal_partial, Fs);


 %%
nfft = 2^nextpow2(N);    % Next power of 2 for efficient FFT
f = Fs*(0:(nfft/2))/nfft; % Frequency vector for plotting

% FFT of noisy speech
D_fft = fft(d, nfft);
D_mag = abs(D_fft/nfft);           % Normalize
D_mag_single = D_mag(1:nfft/2+1);  % Single-sided spectrum
D_mag_single(2:end-1) = 2*D_mag_single(2:end-1); % Scale except DC and Nyquist

% FFT of output signal
Out_fft = fft(output_signal_partial, nfft);
Out_mag = abs(Out_fft/nfft);
Out_mag_single = Out_mag(1:nfft/2+1);
Out_mag_single(2:end-1) = 2*Out_mag_single(2:end-1);

% Plot
figure;
subplot(2,1,1);
plot(f, D_mag_single, 'b', 'LineWidth', 1.2);
title('FFT of Noisy Speech');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

subplot(2,1,2);
plot(f, Out_mag_single, 'r', 'LineWidth', 1.2);
title('FFT of Output Signal After Noise Suppression');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
