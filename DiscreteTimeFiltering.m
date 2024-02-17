% ELE409 Preliminary Work 6
% Mehmet Nurettin Gündüz - 21990887

% Q1
% 
% function y= dftfilt(x,h,N)
% X = fft(x,N);
% H = fft(h,N);
% Y = X .* H;
% y = ifft(Y,N);
% end

% Q2
% 
% function y = convfilt(x,h,N)
% x = [x  zeros(1, N-length(x))];
% h = [h  zeros(1, N-length(h))];
% X = fft(x,N);
% H = fft(h,N);
% Y = X .* H;
% y = ifft(Y,N);
% end

% Q3
% n = 0 : 255;
% x1 = cos(0.1*pi*n) + cos(0.3*pi*n) + cos(0.5*pi*n) + cos(0.7*pi*n) + cos(0.9*pi*n);

% Q4
% 
% h1 = fir1(50, 0.2, 'low');
% y = conv(x1, h1);
% y1 = dftfilt(x1, h1, 256);
% y2 = convfilt(x1, h1, 256); 
% y3 = dftfilt(x1, h1, 512); 
% y4 = convfilt(x1, h1, 512);

% Q5
% 
% subplot(5,1,1)
% plot(y);
% title('y');
%  
% subplot(5,1,2)
% plot(y1);
% title('y1');
%  
% subplot(5,1,3)
% plot(y2);
% title('y2');
%  
% subplot(5,1,4)
% plot(y3);
% title('y3');
%  
% subplot(5,1,5)
% plot(y4);
% title('y4');
% 
% [Y,wy] = freqz(y);
% [Y1,wy1] = freqz(y1);
% [Y2,wy2] = freqz(y2);
% [Y3,wy3] = freqz(y3);
% [Y4,wy4] = freqz(y4);
%  
% subplot(5,1,1);
% plot(wy/pi,abs(Y))
% title('Y');
% 
% subplot(5,1,2);
% plot(wy1/pi,abs(Y1))
% title('Y1');
% 
% subplot(5,1,3);
% plot(wy2/pi,abs(Y2))
% title('Y2');
% 
% subplot(5,1,4);
% plot(wy3/pi,abs(Y3))
% title('Y3');
% 
% subplot(5,1,5);
% plot(wy4/pi,abs(Y4))
% title('Y4');

% Q6
% h1 = fir1(50, 0.2, 'low');
%[x,fs] = audioread('sound.wav');
% y = dftfilt(transpose(x), h1, length(x) + length(h1) - 1);
% [h,w] = freqz(y);
% 
% subplot(2,1,1);
% plot(y)
% title('Output Signal');
% 
% subplot(2,1,2);
% plot(w/pi, abs(h));
% title('Spectrum of Output Signal');

% Q1
% 
% function y= dftfilt(x,h,N)
% X = fft(x,N);
% H = fft(h,N);
% Y = X .* H;
% y = ifft(Y,N);
% end

Q2
% 
% function y = convfilt(x,h,N)            
% x = [x  zeros(1, N-length(x))];
% h = [h  zeros(1, N-length(h))];
% X = fft(x,N);
% H = fft(h,N);
% Y = X .* H;
% y = ifft(Y,N);
% end

