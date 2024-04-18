% Define filter coefficients for the input signal
b = [1 -0.4944 0.64];
a = [1 -1.3335 0.49];
x = [1 1 1 1];  % Define the input signal
L = 20;  % Length of the output signal

% Call the inout function to compute the output signal
y = inout(b, a, x, L);

% Function to compute the output signal using a difference equation
function y = inout(b, a, x, L)
    % Initialize the output signal y
    y = zeros(1, L);
    % Loop over each sample index n up to length L
    for n = 1 : L
        % Initialize the sums for input and output terms
        Sum_X = 0;
        Sum_Y = 0;
        % Loop over the coefficients a(k) except a(1)
        for N = 2 : length(a)
            % Check if the index is within the signal range
            if (n-N+1) > 0
                % Update the output sum with a(k)*y(n-k+1)
                Sum_Y = Sum_Y + a(N) * y(n-N+1);
            end
        end
        % Loop over the coefficients b(k)
        for M = 1 : length(b)
            % Check if the index is within the signal and filter ranges
            if (n-M+1) > 0 && (n-M+1) < length(x)
                % Update the input sum with b(k)*x(n-k+1)
                Sum_X = Sum_X + b(M) * x(n-M+1);
            end
        end
        % Compute the output sample y(n) using the difference equation
        y(n) = (Sum_X - Sum_Y) / a(1);
    end
end

% Define the myAR function to recursively implement an autoregressive filter
function y = myAR(x, a)
    N = 1 : length(a);
    if length(x) > 1
        % Compute the output sample y using the recursive formula
        y = x - a(N) * myAR(x(length(x)-1), a-N);
    end
end

% Generate a Gaussian random signal of length 1024
x = randn(1, 1024);
% Compute the magnitude spectrum of the Gaussian signal using FFT
Magnitude = fftshift(abs(fft(x)));
f = linspace(-1, 1, length(Magnitude));

% Plot the Gaussian signal and its magnitude spectrum
subplot(2,1,1)
plot(x);
title('Gaussian signal');
subplot(2,1,2)
plot(f, Magnitude);
title('Magnitude spectrum of Gaussian signal');

% Define filter coefficients for two different filter designs
b1 = [1 -0.4944 0.64];
a1 = [1 -1.3335 0.49];
b2 = [1 0.4944 0.64];
a2 = [1 1.3335 0.49];

% Apply the filters to the Gaussian signal and compute magnitude spectra
y_a1 = filter(b1, a1, x);
y_a2 = filter(b2, a2, x);
Magnitude2 = fftshift(abs(fft(y_a1)));
Magnitude3 = fftshift(abs(fft(y_a2)));

% Plot the original signal, output from filter h1, and output from filter h2
subplot(3,1,1)
plot(Magnitude)
title("Gaussian signal")
subplot(3,1,2)
plot(Magnitude2)
title("Output from h1")
subplot(3,1,3)
plot(Magnitude3)
title("Output from h2")

% Apply the filters to the Gaussian signal and sum the outputs
y_b1 = filter(b1, a1, x);
y_b2 = filter(b2, a2, x);
y_b3 = y_b1 + y_b2;
Magnitude1 = fftshift(abs(fft(y_b3)));

% Plot the original signal and the output from parallel connection
subplot(2,1,1)
plot(Magnitude)
title("Gaussian signal")
subplot(2,1,2)
plot(Magnitude1)
title("Output from parallel connection")

% Apply the filters in cascade connection to the Gaussian signal
y_c1 = filter(b1, a1, x);
y_c2 = filter(b2, a2, y_c1);
Magnitude = fftshift(abs(fft(x)));
Magnitude1 = fftshift(abs(fft(y_c2)));

% Plot the original signal and the output from cascade connection
subplot(2,1,1)
plot(Magnitude)
title("Gaussian signal")
subplot(2,1,2)
plot(Magnitude1)
title("Output from cascade connection")

% Define complex exponential sequences for different frequencies and compute frequency responses
e_a = [0.8 * exp(1i * 0.1 * pi) 0.8 * exp(-1i * 0.1 * pi)];
freqz(1, poly(e_a), 100);
e_b = [0.8 * exp(1i * 0.5 * pi) 0.8 * exp(-1i * 0.5 * pi)];
freqz(1, poly(e_b), 100);
e_c = [0.8 * exp(1i * 0.9 * pi) 0.8 * exp(-1i * 0.9 * pi)];
freqz(1, poly(e_c), 100);
e_d =[0.1* exp(1i * 0.5 * pi) 0.1 * exp(-1i * 0.5 * pi)];
freqz(1, poly(e_d), 100);
e_e = [0.9 * exp(1i * 0.5 * pi)  0.9 * exp(-1i * 0.5 * pi)];
freqz(1, poly(e_e), 100);

% Define filter coefficients and plot impulse response, step response, and pole-zero diagram
pay_a = [1 0.7264 0.64];
payda_a = [1 -0.6356 0.49];
subplot(3,1,1);
impz(pay_a, payda_a, 30);
title("Impulse Response")
subplot(3,1,2);
stepz(pay_a, payda_a, 30);
title("Step Response")
subplot(3,1,3);
zplane(pay_a, payda_a);
title("Pole-Zero Diagram")

% Repeat for other filter coefficients
pay_b = [1 1.1350 1.5625];
payda_b = [1 -0.6356 0.49];
subplot(3,1,1);
impz(pay_b, payda_b, 30);
title("Impulse Response")
subplot(3,1,2);
stepz(pay_b, payda_b, 30);
title("Step Response")
subplot(3,1,3);
zplane(pay_b, payda_b);
title("Pole-Zero Diagram")

pay_c = [1 0.7264 0.64];
payda_c = [1 -1.362 2.25];
subplot(3,1,1);
impz(pay_c, payda_c, 30);
title("Impulse Response")
subplot(3,1,2);
stepz(pay_c, payda_c, 30);
title("Step Response")
subplot(3,1,3);
zplane(pay_c, payda_c);
title("Pole-Zero Diagram")

% Define filter coefficients and plot magnitude responses

% Filter coefficients for filter A
pay_a = [1 0.7264 0.64];
payda_a = [1 -0.6356 0.49];
subplot(2,1,1)
freqz(pay_a, payda_a, 100); % Compute and plot the magnitude response
title('Magnitude response for filter A');

% Filter coefficients for filter B
pay_b = [1 1.135 1.5625];
payda_b = [1 -0.6356 0.49];
subplot(2,1,2)
freqz(pay_b, payda_b, 100); % Compute and plot the magnitude response
title('Magnitude response for filter B');

