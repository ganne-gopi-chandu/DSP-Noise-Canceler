function snr_value = compute_snr(original_signal, noisy_signal)
    noise = noisy_signal - original_signal;
    
    signal_power = mean(original_signal.^2);
    
    noise_power = mean(noise.^2);
    
    snr_value = 10 * log10(signal_power / noise_power);
    
    fprintf('SNR: %.2f dB\n', snr_value);
end

s=load("clean_speech.txt");
w=load("external_noise.txt");
d=load("noisy_speech.txt");
Fs=44100;

snr_value= compute_snr(s,d);

N = length(d); 
M=round((N));
L = 2; 
epsilon = 1e-6; 
mu = 0.004;         
K=3;
h = zeros(L,1);
estimated_noise = zeros(N,1); 
error_signal = zeros(N,1);    
output_signal = zeros(N,1);   
buffer = zeros(L,1);
for  K =30: 30
   %h = zeros(L,1);
for n = 1:M
    buffer = [w(n); buffer(1:end-1)];
       %h = (buffer/(epsilon + norm(buffer)^2)) * d(n);
    for k = 1:K
        noise_est = h' * buffer;
        e_k = d(n) - noise_est;
        norm_input = norm(buffer)^2;
        h = h + (mu / (epsilon + norm_input)) * e_k * buffer;
        disp(e_k)
    end
   
    estimated_noise(n) = h' * buffer;
    error_signal(n) = d(n) - estimated_noise(n);
    output_signal(n) = error_signal(n);
end

snr_value_final=compute_snr(s(1:M),output_signal(1:M));
end
%sound(output_signal,Fs);
