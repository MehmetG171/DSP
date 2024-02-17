% ELE 409 Preliminary Work 3
% Meehmet Nurettin Gunduz - 21990887

%Question 1.a

% x = SinSamples(3,2*pi,0,2,2*pi*50);
% 
% y= Midrise_Quantizer(x,3); 
% 
% error = x-y; 
% 
% subplot(3,1,1)
% plot(x)
% title('Original Signal') 
% 
% subplot(3,1,2)
% plot(y)
% title('Quantized Version') 
% 
% subplot(3,1,3)
% plot(error)
% title('Quantization Error') 
% 
% SNR = 10 * log10(var(x) / var(error)); 

% %Quesion 1.b
% 
% x = SinSamples(3,2*pi,0,2,2*pi*50); 
% 
% y = Midtread_Quantizer(x,3); 
% 
% error = x-y; 
% 
% subplot(3,1,1)
% plot(x)
% title('Original Signal') 
% 
% subplot(3,1,2)
% plot(y)
% title('Quantized Version')
% 
% subplot(3,1,3)
% plot(error)
% title('Quantization Error') 
% 
% SNR = 10 * log10(var(x) / var(error)); 

% %Question 2
% 
% x = SinSamples(3,2*pi,0,2,2*pi*50);
% 
% x1 = 0.5 * x;
% 
% y1 = Midrise_Quantizer(x1,3); 
% 
% error1 = x1-y1; 
% 
% SNR1 = 10 * log10(var(x1) / var(error1)) 
% 
% x2 = 2 * x; 
% 
% y2 = Midrise_Quantizer(x2,3);  
% 
% error2 = x2 - y2; 
% 
% SNR2 = 10 * log10(var(x2) / var(error2)) 

% %Question 3
% 
% x = audioread('sound3.wav'); 
% 
% y = Midrise_Quantizer(x,3); 
% 
% sound(y) 
% 
% error = x'- y; 
% 
% plot(error) 
% title('Error Signal') 
% 
% SNR = 10 * log10(var(x) / var(error));
% 
% %Question 4

% x = audioread('sound3.wav'); 
% 
% y = Midrise_Quantizer(x,4); 
% 
% sound(y) 
% 
% error = x'-y; 
% 
% plot(error) 
% title('Error Signal') 
% 
% SNR = 10 * log10(var(x) / var(error)) 
% 
% function x = SinSamples(A,w,teta,d,ws)
% T = 0 : (2*pi) / ws : d;
% x = A * sin(w*T + teta);
% % stem(T,x)
% end
% 
% function Q = Midrise_Quantizer(x,bit)
% 
% x_m = max(abs(x));
% 
% delta = 2*x_m / (2^bit);
% 
%  for n = 1:length(x)
% 
%     k = abs(floor(x(n) / delta)); 
%     
%     if k < 2^(bit-1)  
%         Q(n) = ( (floor(x(n) / delta)) + 0.5 ) * delta;
%     
%     elseif x(n) > 0       
%         Q(n) = ( (2^(bit-1) - 1)  + 0.5 ) * delta;
%     
%     elseif x(n) < 0       
%         Q(n) = ( (-(2^(bit-1)) + 1) - 0.5 ) * delta;
%     end   
%   end
% end
% 
% function Q = Midtread_Quantizer(x,bit)

x_m = max(abs(x));

delta =  2*x_m / (2^bit);
 
for n= 1: length(x)

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

