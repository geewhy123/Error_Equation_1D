function taunew = tefft(x,tau,F)

% close all
% clear all
Fs = 2*(length(x)-2)%1000;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = 1000;                     % Length of signal
% x = (0:L-1)*T;



% y = sin(2*pi*x)+0.1*sin(100*pi*x);

tau0 = tau;
% load('tau.mat')
tau = F;
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


% % butterworth filter?
% % % t = 1000;
% % % a = T/t;
% % % yf = filter(a,[1 a-1],y);
% % butter(y,0.1)
% % figure
% % plot(x,yf)
% % % ,x,tau,'*')
% b = 1;
% a = [1 -1];
% % s3_f = filter(b,a,y);
% filtLow = 60;64;
% Fs = N;1024;
% [b, a] = butter(10, filtLow/(Fs/2), 'low');
% s3_f = filtfilt(b, a, y);
% % s3_f = filtfilt([0.2 0.2 0.2 0.2 0.2],1,y)
% figure
% plot(x,s3_f)
% s3_f_fft = abs(fft(s3_f));
% s3_f_fft = s3_f_fft(1:end/2+1);
% figure
% plot(f,s3_f_fft)
% taunew(2:N-1) = s3_f;
% taunew(N) = 0;
% % taunew = tau0;
% % save('tefilt.mat','taunew')

%discard?
Fs = 2*(N-2);%round(1/diff(x(2:3)))
Fs = N-2;
Y = fft(tau);%tau0);
L = length(Y);
% s3_f_fft(end/10:end) = 0;
f = Fs/2*linspace(0,1,NFFT/2+1)%
f = Fs*linspace(0,1,L)
cutoff = 100;
NY = Y;
% NY(f>cutoff & f<(Fs-cutoff))=0;
NY(f>cutoff )=0;
Ny = real(ifft(NY));
% taunew = ifft(s3_f_fft);
% taunew = real(taunew);
figure
% plot(f,s3_f_fft(1:length(f)));
length(f)
length(NY(1:NFFT/2+1))
% plot(f,2*abs(NY(1:NFFT/2+1)))%abs(NY))
figure
plot(x,Ny(2:end-1))
taunew = (Ny);
% error('1')

%subtract mean
% taunew = tau0-mean(tau0);
% figure
% length(x)
% length(taunew)
% plot(x,taunew(2:length(x)+1))
end