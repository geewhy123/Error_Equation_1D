function taunew = tefft(x,tau,F)

% taunew = F;
% return;

% close all
% clear all
Fs = 1*(length(x)-2)%1000;                    % Sampling frequency
% T = 1/Fs;                     % Sample time
% L = 1000;                     % Length of signal
% x = (0:L-1)*T;

% F = ones(size(F))
% F = sin(2*pi*x)+0.3*sin(40*pi*x)+0.1*cos(100*pi*x);
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
% Y = fft(y)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
% f = Fs/2*linspace(0,1,L);
% Plot single-sided amplitude spectrum.
figure
subplot(2,2,1)
plot(x,y,'*-');
xlabel('x')
ylabel('residual source(2,q,4)')
% figure
subplot(2,2,2)
stem(f,2*abs(Y(1:NFFT/2+1))) 
MM = 2*max(abs(Y));
ylim([0 MM])
% plot(f,abs(Y(2:end/2)));
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

% error('1')
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
% filtLow = 64;
% Fs = N;1024;
% % [b, a] = butter(10, filtLow/(Fs/2), 'low');
% [b, a] = butter(2, filtLow/(Fs/2));%, 'low');
% % s3_f = filtfilt(b, a, y);
% s3_f = filter(b,a,y);
% % s3_f = filtfilt([0.2 0.2 0.2 0.2 0.2],1,y)
% figure
% plot(x,s3_f)
% s3_f_fft = abs(fft(s3_f));
% s3_f_fft = s3_f_fft(1:end/2+1);
% figure
% plot(f,s3_f_fft)
% taunew(2:N-1) = s3_f;
% taunew = [NaN;taunew';NaN]';
% % taunew = tau0;
% % % save('tefilt.mat','taunew')
% return;


% %gaussian
% fn = y;
% m = 101;
% mu = 6;
% n = L;
% h = build_gaussian_filter(m,mu/(4*n),n);
% fnh = perform_convolution(fn,h);
% subplot(223);plot(x,fnh,'-*');
% YY = fft(fnh,NFFT)/L;
% subplot(224);stem(f,2*abs(YY(1:NFFT/2+1))) ;
% ylim([0 MM])
% taunew = [NaN;fnh;NaN];
% return;






% error('1')
% x = [0; x ;0];
% tau = sin(2*pi*x)+sin(40*pi*x);
% x = x(2:end-1);
tau = tau(2:end-1);

%discard?
Fs = 1*(N-2);%round(1/diff(x(2:3)))
% Fs = N-2;
Y = fft(tau);%tau0);
L = length(Y);
% s3_f_fft(end/10:end) = 0;
f = Fs/2*linspace(0,1,NFFT/2+1)%
% f = Fs*linspace(0,1,L)
cutoff = 200;
NY = Y;
% NY(f>cutoff & f<(Fs-cutoff))=0;
% NY(f>cutoff & f<Fs/4)=0;
% NY(f<(Fs/2-cutoff) & f>Fs/4)=0;

 NY(cutoff+1:end-cutoff+1) = 0;
%  phi = 0.1;
%  NY(2) = NY(2)*exp(j*phi);
%  NY(end)=NY(end)*exp(-j*phi);
% NY(1) = 0;

% NY
% error('1')
% % NY(1:cutoff+1)=0;
% % NY(end-cutoff-1:end) = 0;

% % band-stop
% cutoffL = 20;
% cutoffR = 58;
% NY(cutoffL+1:cutoffR) = 0;
% NY(end-cutoffR-1:end-cutoffL) = 0;


% f> cutoff
% f > 3*Fs/4
% f< (Fs/2-cutoff)
% NY

% NY(f>cutoff )=0;
% NY(f<(Fs-cutoff))=0;
% NY
Ny = real(ifft(NY));
% error('1')
% taunew = ifft(s3_f_fft);
% taunew = real(taunew);
% figure
% plot(f,s3_f_fft(1:length(f)));
length(f)
length(NY(1:NFFT/2+1))
2*abs(NY(1:NFFT/2+1))
subplot(2,2,4)
stem(f,2*abs(NY(1:NFFT/2+1))/L)%abs(NY))
ylim([0 MM])
% figure
% error('1')
subplot(2,2,3)
plot(x,Ny,'*-')%(2:end-1))
taunew = [NaN;Ny;NaN];

% mean(Ny)
% error('1')
% taunew = taunew+0.002;
% error('1')

%subtract mean
% taunew = tau-mean(tau);
% figure
% length(x)
% length(taunew)
% plot(x,taunew)
% % plot(x,taunew(2:length(x)+1))
% taunew = [NaN; taunew ; NaN];
end