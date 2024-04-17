
% Generate sinusoidal samples
x1 = SinSamples(1, 2*pi*1000, pi/6, 2e-3, 2*pi*16000); % Sinusoid with phase shift
X1 = Spectrum(x1, 512, 2*pi*16000); % Compute and plot spectrum

x2 = SinSamples(1, 2*pi*1000, 0, 2e-3, 2*pi*16000); % Sinusoid without phase shift
X2 = Spectrum(x2, 512, 2*pi*16000);

x3 = SinSamples(1, 2*pi*1000, 0, 10e-3, 2*pi*4000); % Lower frequency sinusoid
X3 = Spectrum(x3, 512, 2*pi*4000);

x4 = SinSamples(1, 2*pi*5000, 0, 10e-3, 2*pi*4000); % Higher frequency sinusoid
X4 = Spectrum(x4, 512, 2*pi*4000);

x5 = SinSamples(1, 2*pi*1000, 0, 1e-3, 2*pi*16000); % Short duration sinusoid
X5 = Spectrum(x5, 512, 2*pi*16000);

x6 = SinSamples(1, 2*pi*1000, 0, 15e-4, 2*pi*16000); % Very short duration sinusoid
X6 = Spectrum(x6, 512, 2*pi*16000);

x7 = SinSamples(1, 2*pi*1000, 0, 2e-3, 2*pi*16000); % Normal duration sinusoid
X7 = Spectrum(x7, 512, 2*pi*16000);

x8 = SinSamples(3, 2*pi*5000, pi/6, 2e-3, 2*pi*16000); % High-amplitude sinusoid

% Combine signals and compute spectrum
xs = x1 + x8; % Add two signals
Xs = Spectrum(xs, 512, 2*pi*16000); % Compute and plot spectrum

xm = x1 .* x8; % Multiply two signals
Xm = Spectrum(xm, 512, 2*pi*16000); % Compute and plot spectrum

% Perform circular shift and compute spectrum
x9 = circshift(x5, -1); % Circular shift
X9 = Spectrum(x9, 512, 2*pi*16000); % Compute and plot spectrum

% Read audio file and compute spectrum
[y, Fs] = audioread('sound1.wav', [1, 512]); % Read audio file
Y = fft(y, 512); % Compute FFT of audio
time = (1 : 1 : 512) / Fs; % Time vector
frequency = (-512/2 : 1 : 512/2 - 1) * Fs / 512; % Frequency vector

% Plot audio signal and its spectrum
plot(time, y) % Plot audio signal over time
plot(frequency, fftshift(abs(Y))) % Plot spectrum

% Generate sinc function and compute spectrum
n = 0 : 1 : 255; % Time index
x10 = sinc(0.2 * (n - 128)); % Sinc function
plot(x10) % Plot sinc function

X10 = fftshift(abs(fft(x10))); % Compute and plot spectrum
plot(X10)

xdown1 = downsample(x10, 2); % Downsample sinc function
plot(xdown1) % Plot downsampled signal

XDOWN1 = fftshift(abs(fft(xdown1))); % Compute and plot spectrum of downsampled signal
plot(XDOWN1)

xinterp1 = interp(xdown1, 2); % Interpolate downsampled signal
plot(xinterp1) % Plot interpolated signal

XINTERP1 = fftshift(abs(fft(xinterp1))); % Compute and plot spectrum of interpolated signal
plot(XINTERP1)

% Generate another sinc function and compute spectrum
n = 0 : 1 : 255; % Time index
x11 = sinc(0.8 * (n - 128)); % Sinc function
plot(x11) % Plot sinc function

X11 = fftshift(abs(fft(x11))); % Compute and plot spectrum
plot(X11)

xdown2 = downsample(x11, 2); % Downsample sinc function
plot(xdown2) % Plot downsampled signal

XDOWN2 = fftshift(abs(fft(xdown2))); % Compute and plot spectrum of downsampled signal
plot(XDOWN2)

xinterp2 = interp(xdown2, 2); % Interpolate downsampled signal
plot(xinterp2) % Plot interpolated signal

XINTERP2 = fftshift(abs(fft(xinterp2))); % Compute and plot spectrum of interpolated signal
plot(XINTERP2)

% Function to generate sinusoidal samples
function x = SinSamples(A, w, teta, d, ws)
    T = 0 : (2 * pi) / ws : d; % Time vector
    x = A * sin(w * T + teta); % Generate sinusoid
end

% Function to compute spectrum using FFT and plot it
function X = Spectrum(x, n, ws)
    X = fft(x, n); % Compute FFT
    k = (-n/2 : 1 : n/2 - 1) * ws / n; % Frequency vector
    plot(k, fftshift(abs(X))) % Plot spectrum
end
