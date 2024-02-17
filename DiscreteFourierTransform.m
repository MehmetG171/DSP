% ELE 409 - DSP LAB - PRELIMINARY WORK 1
% MEHMET NURETTIN GUNDUZ - 21990887

%  x1 = SinSamples(1, 2*pi*1000, pi/6, 2e-3, 2*pi*16000);
% 
% X1 = Spectrum(x1,512,2*pi*16000);

% ------

% x2 = SinSamples(1, 2*pi*1000, 0, 2e-3, 2*pi*16000);
% 
% X2 = Spectrum(x2,512,2*pi*16000);
% 
% ------

% x3 = SinSamples(1, 2*pi*1000, 0, 10e-3, 2*pi*4000);
% 
% X3 = Spectrum(x3,512,2*pi*4000);

% ------

% x4 = SinSamples(1, 2*pi*5000, 0, 10e-3, 2*pi*4000);
% 
% X4 = Spectrum(x4,512,2*pi*4000);

% ------

% x5 = SinSamples(1, 2*pi*1000, 0, 1e-3, 2*pi*16000);
% 
% X5 = Spectrum(x5,512,2*pi*16000);

% ------

% x6 = SinSamples(1, 2*pi*1000, 0, 15e-4, 2*pi*16000);
% 
% X6 = Spectrum(x6,512,2*pi*16000);

% ------

% x7 = SinSamples(1, 2*pi*1000, 0, 2e-3, 2*pi*16000);
% 
% X7 = Spectrum(x7,512,2*pi*16000);
% 
% ------

%  x8 = SinSamples(3, 2*pi*5000, pi/6, 2e-3, 2*pi*16000);
% 
%  xs = x1 + x8;
%  
%  Xs = Spectrum(xs,512,2*pi*16000);
% 
%  xm = x1 .* x8;
% 
%  Xm = Spectrum(xm,512,2*pi*16000);
% 
% ------

% x9 = circshift(x5,-1);
%  
% X9 = Spectrum(x9,512,2*pi*16000);

% ------

% [y,Fs] = audioread('sound1.wav',[1,512]);
% Y = fft(y,512);
% time = (1 : 1 : 512) / Fs;
% frequency = (-512/2 : 1 : 512/2 - 1) * Fs / 512;

% plot(time,y)
% plot(frequency, fftshift(abs(Y)))

% ------

% n = 0 : 1 : 255;
% x10 = sinc(0.2 * (n-128)); 
% plot(x10);

% X10 = fftshift (abs (fft (x10)));
% plot(X10);

% xdown1 = downsample (x10,2); 
% plot(xdown1);

% XDOWN1 = fftshift (abs (fft (xdown1))); 
% plot(XDOWN1);
% 
% xinterp1 = interp (xdown1,2); 
% plot(xinterp1);
% 
% XINTERP1 = fftshift (abs (fft (xinterp1))); 
% plot(XINTERP1);

% ------

% n = 0 : 1 : 255;
% x11 = sinc(0.8 * (n-128));
% plot(x11);
% 
% X11 = fftshift (abs (fft (x11)));
% plot(X11);
% 
% xdown2 = downsample (x11,2); 
% plot(xdown2);
% 
% XDOWN2 = fftshift (abs (fft (xdown2))); 
% plot(XDOWN2);

% xinterp2 = interp (xdown2,2);
% plot(xinterp2);
% 
% XINTERP2 = fftshift (abs (fft (xinterp2))); 
% plot(XINTERP2);

% ------

function x = SinSamples(A,w,teta,d,ws)
T = 0 : (2*pi) / ws : d;
x = A * sin(w*T + teta);
% stem(T,x)
end

function X = Spectrum(x,n,ws)
X = fft(x,n);
k = (-n/2 : 1 : n/2 - 1) * ws / n;
plot(k,fftshift(abs(X)))
end