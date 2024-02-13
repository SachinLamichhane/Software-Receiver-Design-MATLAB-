%Calculating power spectral density

% Sample MATLAB code for calculating Power Spectral Density (PSD)

% Step 1: Generate or load time-domain signal

t = 0:1/f_s:1; % Time vector
signal = filter_signal;
% Step 2: Apply FFT
N = length(signal);
fft_result = fft(signal);
%figure('Name','Periodogram'), periodogram(demod_sig);
%figure('Name','Periodogram match'), periodogram(filter_signal);
% Step 3: Calculate PSD
psd = (1/(f_s*N)) * abs(fft_result).^2;

% Step 4: Frequency axis
f = f_s*(0:(N/2))/N;

% Step 5: Plotting
figure;
plot(f, 10*log10(psd(1:N/2+1)));
title('Power Spectral Density');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;