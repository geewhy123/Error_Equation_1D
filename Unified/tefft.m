close all
clear all
Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = 1000;                     % Length of signal
x = (0:L-1)*T;



% y = sin(2*pi*x)+0.1*sin(100*pi*x);


load('tau.mat')
N = length(x);
x= x(2:end-1);

% tau = sin(2*pi*x);
L = length(x);
y = tau;
% y = sin(2*pi*x)+0.1*sin(400*pi*x);
% y = rand(size(x))-0.5;
% y = randn(size(x));
% y = x;
% y(1:ceil(end/2)) = 1;
% y(ceil(end/2)+1:end) = -1;
y = y(2:end-1);

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure
plot(x,y);
figure
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

% % t = 1000;
% % a = T/t;
% % yf = filter(a,[1 a-1],y);
% butter(y,0.1)
% figure
% plot(x,yf)
% % ,x,tau,'*')

b = 1;
a = [1 -1];
% s3_f = filter(b,a,y);
filtLow = 10;
Fs = N;1024;
[b, a] = butter(10, filtLow/(Fs/2), 'low');
s3_f = filtfilt(b, a, y);
% s3_f = filtfilt([0.2 0.2 0.2 0.2 0.2],1,y)
figure
plot(x,s3_f)
s3_f_fft = abs(fft(s3_f));
s3_f_fft = s3_f_fft(1:end/2+1);
figure
plot(f,s3_f_fft)
taunew(2:N-1) = s3_f;
taunew(N) = 0;
save('tefilt.mat','taunew')