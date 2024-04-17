% Generate a sinusoidal signal
x = SinSamples(3, 2*pi, 0, 2, 2*pi*50);

% Quantize the signal using midrise quantization with 3 bits
y_midrise = Midrise_Quantizer(x, 3); 

% Calculate the quantization error
error_midrise = x - y_midrise; 

% Plot the original signal, quantized version, and quantization error
subplot(3, 1, 1)
plot(x)
title('Original Signal') 

subplot(3, 1, 2)
plot(y_midrise)
title('Quantized Version (Midrise)') 

subplot(3, 1, 3)
plot(error_midrise)
title('Quantization Error (Midrise)') 

% Calculate the Signal-to-Noise Ratio (SNR) in dB
SNR_midrise = 10 * log10(var(x) / var(error_midrise)); 

% Generate the same sinusoidal signal
x = SinSamples(3, 2*pi, 0, 2, 2*pi*50); 

% Quantize the signal using midtread quantization with 3 bits
y_midtread = Midtread_Quantizer(x, 3); 

% Calculate the quantization error
error_midtread = x - y_midtread; 

% Plot the original signal, quantized version, and quantization error
subplot(3, 1, 1)
plot(x)
title('Original Signal') 

subplot(3, 1, 2)
plot(y_midtread)
title('Quantized Version (Midtread)') 

subplot(3, 1, 3)
plot(error_midtread)
title('Quantization Error (Midtread)') 

% Calculate the Signal-to-Noise Ratio (SNR) in dB
SNR_midtread = 10 * log10(var(x) / var(error_midtread)); 

% Generate a sinusoidal signal with reduced amplitude
x1 = 0.5 * x;

% Quantize the reduced amplitude signal using midrise quantization with 3 bits
y1_midrise = Midrise_Quantizer(x1, 3); 

% Calculate the quantization error
error1_midrise = x1 - y1_midrise; 

% Calculate the Signal-to-Noise Ratio (SNR) for the reduced amplitude signal
SNR1_midrise = 10 * log10(var(x1) / var(error1_midrise)); 

% Generate a sinusoidal signal with increased amplitude
x2 = 2 * x; 

% Quantize the increased amplitude signal using midrise quantization with 3 bits
y2_midrise = Midrise_Quantizer(x2, 3); 

% Calculate the quantization error
error2_midrise = x2 - y2_midrise; 

% Calculate the Signal-to-Noise Ratio (SNR) for the increased amplitude signal
SNR2_midrise = 10 * log10(var(x2) / var(error2_midrise)); 

% Load an audio signal
x_audio = audioread('sound3.wav'); 

% Quantize the audio signal using midrise quantization with 3 bits
y_audio_midrise = Midrise_Quantizer(x_audio, 3); 

% Play the quantized audio signal
sound(y_audio_midrise) 

% Calculate the quantization error for the audio signal
error_audio_midrise = x_audio' - y_audio_midrise; 

% Plot the quantization error for the audio signal
plot(error_audio_midrise) 
title('Error Signal (Audio, Midrise Quantization)') 

% Calculate the Signal-to-Noise Ratio (SNR) for the audio signal
SNR_audio_midrise = 10 * log10(var(x_audio) / var(error_audio_midrise)); 

% Quantize the audio signal using midrise quantization with 4 bits
y_audio_midrise_4bits = Midrise_Quantizer(x_audio, 4); 

% Play the quantized audio signal
sound(y_audio_midrise_4bits) 

% Calculate the quantization error for the audio signal
error_audio_midrise_4bits = x_audio' - y_audio_midrise_4bits; 

% Plot the quantization error for the audio signal
plot(error_audio_midrise_4bits) 
title('Error Signal (Audio, Midrise Quantization, 4 bits)') 

% Calculate the Signal-to-Noise Ratio (SNR) for the audio signal
SNR_audio_midrise_4bits = 10 * log10(var(x_audio) / var(error_audio_midrise_4bits)); 

% Function to generate a sinusoidal signal
function x = SinSamples(A, w, teta, d, ws)
    T = 0 : (2*pi) / ws : d;
    x = A * sin(w * T + teta);
end

% Function for midrise quantization
function Q = Midrise_Quantizer(x, bit)
    x_m = max(abs(x));
    delta = 2 * x_m / (2^bit);
    Q = zeros(size(x));
    for n = 1:length(x)
        k = abs(floor(x(n) / delta)); 
        if k < 2^(bit-1)  
            Q(n) = ( (floor(x(n) / delta)) + 0.5 ) * delta;
        elseif x(n) > 0       
            Q(n) = ( (2^(bit-1) - 1)  + 0.5 ) * delta;
        elseif x(n) < 0       
            Q(n) = ( (-(2^(bit-1)) + 1) - 0.5 ) * delta;
        end   
    end
end

% Function for midtread quantization
function Q = Midtread_Quantizer(x, bit)
    x_m = max(abs(x));
    delta = 2 * x_m / (2^bit);
    Q = zeros(size(x));
    for n=1:length(x)
        k = abs(round(x(n) / delta)); 
        if k < 2^(bit-1)
            Q(n) = round(x(n) / delta) * delta;
        elseif x(n) > 0       
            Q(n) = (2^(bit-1)-1) * delta;
        elseif x(n) < 0       
            Q(n) = (-(2^(bit-1))) * delta;        
        end    
    end
end
