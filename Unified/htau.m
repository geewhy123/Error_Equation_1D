close all
clear all
load('tauNU.mat')

tau3264 = tau32;
N = 32;
for j = 2:N+1
   h12 = h64(2*(j-1))
   h23 = h64(2*(j-1)+1);
   tau3264(j) = (h12*tau64(2*(j-1))+h23*tau64(2*(j-1)+1))/(h12+h23); 
end
plot(x32,tau32,'-+',x64,tau64,'-*',x32,tau3264,'o')
set(gca,'FontSize',15)
legend('N=32','N=64','N=64 interpolated back to N=32 by weighted average')
xlabel('x')
ylabel('\tau')