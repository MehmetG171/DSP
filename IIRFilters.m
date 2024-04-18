% Load the audio file 'sound4.wav' and plot its waveform
[x_input, fs] = audioread('sound4.wav');
plot(x_input);
title('Waveform');

% Calculate the FFT of the input signal and plot its magnitude response
X_INPUT = fftshift(abs(fft(x_input)));
freq = (-fs/2):(fs/length(x_input)):(fs/2)-(fs/length(x_input));
plot(freq, X_INPUT);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Response');

% Apply Butterworth filter and plot the filtered signal and its magnitude response
x_butterworth = filter(butterworth, x_input);
X_BUTTERWORTH = fftshift(abs(fft(x_butterworth)));
subplot(1,2,1);
plot(x_butterworth);
title('Butterworth');
subplot(1,2,2);
plot(freq, X_BUTTERWORTH);
title('Magnitude Response of Butterworth');

% Apply Chebyshev Type I filter and plot the filtered signal and its magnitude response
y_chebyshev1 = filter(chebyshev1, x_input);
Y_CHEBYSHEV1 = fftshift(abs(fft(y_chebyshev1)));
subplot(1,2,1);
plot(y_chebyshev1);
title('Chebyshev Type I');
subplot(1,2,2);
plot(freq, Y_CHEBYSHEV1);
title('Magnitude Response of Chebyshev Type I');

% Apply Chebyshev Type II filter and plot the filtered signal and its magnitude response
y_chebyshev2 = filter(chebyshev2, x_input);
Y_CHEBYSHEV2 = fftshift(abs(fft(y_chebyshev2)));
subplot(1,2,1);
plot(y_chebyshev2);
title('Chebyshev Type II');
subplot(1,2,2);
plot(freq, Y_CHEBYSHEV2);
title('Magnitude Response of Chebyshev Type II');

% Apply Butterworth bandpass filter and plot the filtered signal and its magnitude response
x_butterworthbp = filter(butterworthbp, x_input);
X_BUTTERWORTHBP = fftshift(abs(fft(x_butterworthbp)));
subplot(1,2,1);
plot(x_butterworthbp);
title('Butterworth Bandpass');
subplot(1,2,2);
plot(freq, X_BUTTERWORTHBP);
title('Magnitude Response of Butterworth Bandpass');

% Apply Chebyshev Type I bandpass filter and plot the filtered signal and its magnitude response
x_chebyshev1bp = filter(chebyshev1bp, x_input);
X_CHEBYSHEV1BP = fftshift(abs(fft(x_chebyshev1bp)));
subplot(1,2,1);
plot(x_chebyshev1bp);
title('Chebyshev Type I Bandpass');
subplot(1,2,2);
plot(freq, X_CHEBYSHEV1BP);
title('Magnitude Response of Chebyshev Type I Bandpass');

% Apply Chebyshev Type II bandpass filter and plot the filtered signal and its magnitude response
x_chebyshev2bp = filter(chebyshev2bp, x_input);
X_CHEBYSHEV2BP = fftshift(abs(fft(x_chebyshev2bp)));
subplot(1,2,1);
plot(x_chebyshev2bp);
title('Chebyshev Type II Bandpass');
subplot(1,2,2);
plot(freq, X_CHEBYSHEV2BP);
title('Magnitude Response of Chebyshev Type II Bandpass');

% Load the audio file 'bana.wav', play it, and filter it using the designed filter Hd
[x_input, fs] = audioread('bana.wav');
sound(x_input);
y = filter(Hd, x_input);
sound(y);

% Convert the second-order sections (SOS) to transfer function (TF) coefficients
[b, a] = sos2tf(SOS);
