% ELE 409 Preliminary Work 5
% Mehmet Nurettin Gunduz - 21990887

% Question 1
%
% rec(40);
% title("Function Test with rec(40)")
% bart(40);
% title("Function Test with bart(40)")
% bla(40);
% title("Function Test with bla(40)")
%
% Functions are in the bottom of the script.

% Question 2
%
% rec_40 = rec(40);
% freqz(rec_40,100);
% title('M = 40');
%
% rec_80 = rec(80);
% freqz(rec_80,100);
% title('M = 80');
%
% bart_40 = bart(40);
% freqz(bart_40,100);
% title('M = 40');
% 
% bart_80 = bart(80);
% freqz(bart_80,100);
% title('M = 80');
% 
% bla_40 = bla(40);
% freqz(bla_40,100);
% title('M = 40');
%
% bla_80 = bla(80);
% freqz(bla_80,100);
% title('M = 80');

% Question 3
%
% for n = 1 : 256
%      x1(n) = cos(2*pi*0.242*n) + cos(2*pi*0.258*n);
%      x2(n) = 2*cos(2*pi*0.25*n) + cos(2*pi*0.008*n);
%      x3(n) = cos(2*pi*0.29*n);
% end
%
% plot(x1);
% title('Time Waveform of x1[n]');
%
% freqz(x1)
% title('Magnitude Response of x1[n]');
% 
% plot(x2);
% title('Time Waveform of x2[n]');
% 
% freqz(x2)
% title('Magnitude Response of x2[n]');

% Question 4
% 
% rec_40_new = [rec_40  zeros(1,(length(x1)-(length(rec_40))))];   
% x1_rec40 = rec_40_new .* x1; 
% freqz(x1_rec40)
% title('x1 through Rectangular Window with M = 40');
% 
% rec_80_new = [rec_80  zeros(1,(length(x1)-(length(rec_80))))];    
% x1_rec80 = rec_80_new .* x1; 
% freqz(x1_rec80)
% title('x1 through Rectangular Window with M = 80');

% Question 5
%
% rec3 = rec(length(x3));
% rec_x3=(x3).*(rec3);
% freqz(rec_x3);
% title('x3 through Rectangular Window');
%
% bart3 = bart(length(x3));
% bart_x3 = (x3) .* (bart3);
% freqz(bart_x3);
% title('x3 through Barlett Window');

function w = rec(M)
for i = 1 : 1 : M
w(i) = 1;                         
end
plot(w);
end

function w = bart(M)
for i = 1 : 1 : M 
w(i) = 1 - abs( (i - M/2) / ( M/2 ) );
end
plot(w)
end

function w = bla(M)
for i = 1 : 1 : M
w(i) = (0.42) - (0.5) * cos( (2*pi*i) / (M-1)) + (0.08) * cos( (4*pi*i) / (M-1));
end
plot(w)
end