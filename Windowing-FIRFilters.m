rec(40); % Plot the rectangular window with length 40
title("Function Test with rec(40)")

bart(40); % Plot the Bartlett window with length 40
title("Function Test with bart(40)")

bla(40); % Plot the Blackman window with length 40
title("Function Test with bla(40)")

% Generate different windows and analyze their frequency responses
rec_40 = rec(40);
freqz(rec_40,100);
title('M = 40');

rec_80 = rec(80);
freqz(rec_80,100);
title('M = 80');

bart_40 = bart(40);
freqz(bart_40,100);
title('M = 40');

bart_80 = bart(80);
freqz(bart_80,100);
title('M = 80');

bla_40 = bla(40);
freqz(bla_40,100);
title('M = 40');

bla_80 = bla(80);
freqz(bla_80,100);
title('M = 80');

% Generate signals x1, x2, x3 and plot their time waveforms
for n = 1 : 256
     x1(n) = cos(2*pi*0.242*n) + cos(2*pi*0.258*n);
     x2(n) = 2*cos(2*pi*0.25*n) + cos(2*pi*0.008*n);
     x3(n) = cos(2*pi*0.29*n);
end

plot(x1);
title('Time Waveform of x1[n]');

plot(x2);
title('Time Waveform of x2[n]');

% Perform frequency analysis on signals x1 and x2
freqz(x1)
title('Magnitude Response of x1[n]');

freqz(x2)
title('Magnitude Response of x2[n]');

% Apply rectangular window to x1 with different lengths (M = 40, 80)
rec_40_new = [rec_40  zeros(1,(length(x1)-(length(rec_40))))];   
x1_rec40 = rec_40_new .* x1; 
freqz(x1_rec40)
title('x1 through Rectangular Window with M = 40');

rec_80_new = [rec_80  zeros(1,(length(x1)-(length(rec_80))))];    
x1_rec80 = rec_80_new .* x1; 
freqz(x1_rec80)
title('x1 through Rectangular Window with M = 80');

% Apply rectangular window to x3
rec3 = rec(length(x3));
rec_x3=(x3).*(rec3);
freqz(rec_x3);
title('x3 through Rectangular Window');

% Apply Bartlett window to x3
bart3 = bart(length(x3));
bart_x3 = (x3) .* (bart3);
freqz(bart_x3);
title('x3 through Bartlett Window');

% Generate a rectangular window function
function w = rec(M)
    for i = 1 : M
        w(i) = 1;  % All values are set to 1, creating a rectangular window
    end
    plot(w);  % Plot the generated window
end

% Generate a Bartlett (triangular) window function
function w = bart(M)
    for i = 1 : M
        w(i) = 1 - abs((i - M/2) / (M/2));  % Calculate the Bartlett window values
    end
    plot(w);  % Plot the generated window
end

% Generate a Blackman window function
function w = bla(M)
    for i = 1 : M
        % Calculate the Blackman window values using the specified formula
        w(i) = 0.42 - 0.5 * cos(2*pi*i / (M-1)) + 0.08 * cos(4*pi*i / (M-1));
    end
    plot(w);  % Plot the generated window
end
