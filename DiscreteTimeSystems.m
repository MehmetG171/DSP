% ELE409 Preliminary Work 2 
% Mehmet Nurettin Gunduz - 21990887

% 2
% b = [1 -0.4944 0.64];
% a = [1 -1.3335 0.49];
% x = [1 1 1 1];
% L = 20;
% y=inout(b,a,x,L);
% 
% function y = inout(b,a,x,L)
%  for n = 1 : L
%     Sum_X = 0;
%     Sum_Y = 0;
%     for N = 2 : length(a)
%         if (n-N+1) > 0
%             Sum_Y = Sum_Y + a(N) * y(n-N+1);
%         end
%     end
%     for M = 1 : length(b)
%        if (n-M+1) > 0
%           if (n-M+1) < length(x)
%             Sum_X = Sum_X + b(M) * x(n-M+1);
%           end
%        end
%     end
%     y(n) = (Sum_X-Sum_Y) / a(1);
%   end
% end

% 3 

% function y = myAR(x,a)
% N = 1 : length(a);
%    if length(x) > 1
%    y = x - a(N) * myAR(x(length(x)-1),a-N);
%    end
% end

% 5
x = randn(1,1024);
% Magnitude = fftshift(abs(fft(x)));
% f = linspace(-1,1,length(Magnitude));
% subplot(2,1,1)
% plot(x);
% title('Gaussian signal');
% subplot(2,1,2)
% plot(f,Magnitude);
% title('Magnitude spectrum of Gaussian signal');

% 6
%  b1 = [1 -0.4944 0.64];
%  a1 = [1 -1.3335 0.49];
%  b2 = [1 0.4944 0.64];
%  a2 = [1 1.3335 0.49];

% 6.a 
%  y_a1 = filter (b1,a1,x);
%  y_a2 = filter (b2,a2,x);
%  Magnitude = fftshift(abs(fft(x)));
%  subplot(3,1,1)
%  plot(Magnitude)
%  title("Gaussian signal")
%  Magnitude2 = fftshift(abs(fft(y_a1)));
%  subplot(3,1,2)
%  plot(Magnitude2)
%  title("Output from h1")
%  Magnitude3 = fftshift(abs(fft(y_a2)));
%  subplot(3,1,3)
%  plot(Magnitude3)
%  title("Output from h2")

% 6.b 
%  y_b1 = filter (b1,a1,x);
%  y_b2 = filter (b2,a2,x);
%  y_b3 = y_b1 + y_b2;
%  Magnitude = fftshift(abs(fft(x)));
%  subplot(2,1,1)
%  plot(Magnitude)
%  title("Gaussian signal")
%  Magnitude1 = fftshift(abs(fft(y_b3)));
%  subplot(2,1,2)
%  plot(Magnitude1)
%  title("Output from parallel connection")

% 6.c 
%  y_c1 = filter (b1,a1,x);
%  y_c2 = filter (b2,a2,y_c1);
%  Magnitude = fftshift(abs(fft(x)));
%  subplot(2,1,1)
%  plot(Magnitude)
%  title("Gaussian signal")
%  Magnitude1 = fftshift(abs(fft(y_c2)));
%  subplot(2,1,2)
%  plot(Magnitude1)
%  title("Output from cascade connection")
 
% 7.a
% e_a = [0.8 * exp(1i * 0.1 * pi) 0.8 * exp(-1i * 0.1 * pi)];
% freqz(1,poly(e_a),100);
% title('7.a');

% 7.b
% e_b = [0.8 * exp(1i * 0.5 * pi) 0.8 * exp(-1i * 0.5 * pi)];
% freqz(1,poly(e_b),100);
% title('7.b');

% 7.c 
% e_c = [0.8 * exp(1i * 0.9 * pi) 0.8 * exp(-1i * 0.9 * pi)];
% freqz(1,poly(e_c),100);
% title('7.c');

% 7.d 
% e_d =[0.1* exp(1i * 0.5 * pi) 0.1 * exp(-1i * 0.5 * pi)];
% freqz(1,poly(e_d),100);
% title('7.d');

% 7.e 
% e_e = [0.9 * exp(1i * 0.5 * pi)  0.9 * exp(-1i * 0.5 * pi)];
% freqz(1,poly(e_e),100);
% title('7.e');

%8.a
% pay_a = [1 0.7264 0.64];
% payda_a = [1 -0.6356 0.49];
% subplot(3,1,1);
% impz(pay_a,payda_a,30);
% title("Impulse Response of a")
% subplot(3,1,2);
% stepz(pay_a,payda_a,30);
% title("Step Response of a")
% subplot(3,1,3);
% zplane(pay_a,payda_a);
% title("Pole-Zero Diagram of a")

%8.b
% pay_b = [1 1.1350 1.5625];
% payda_b = [1 -0.6356 0.49];
% subplot(3,1,1);
% impz(pay_b,payda_b,30);
% title("Impulse Response of b")
% subplot(3,1,2);
% stepz(pay_b,payda_b,30);
% title("Step Response of b")
% subplot(3,1,3);
% zplane(pay_b,payda_b);
% title("Pole-Zero Diagram of b")

% %8.c
% pay_c = [1 0.7264 0.64];
% payda_c = [1 -1.362 2.25];
% subplot(3,1,1);
% impz(pay_c,payda_c,30);
% title("Impulse Response of c")
% subplot(3,1,2);
% stepz(pay_c,payda_c,30);
% title("Step Response of c")
% subplot(3,1,3);
% zplane(pay_c,payda_c);
% title("Pole-Zero Diagram of c")

%9.a
% pay_a = [1 0.7264 0.64];
% payda_a = [1 -0.6356 0.49];
% subplot(2,1,1)
% freqz(pay_a,payda_a,100);
% title('Magnitude response of 8.a');

%9.b
% pay_b = [1 1.135 1.5625];
% payda_b = [1 -0.6356 0.49];
% subplot(2,1,2)
% freqz(pay_b,payda_b,100);
% title('Magnitude response of 8.b');
