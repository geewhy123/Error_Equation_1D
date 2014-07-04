clear all
load('err10.mat');
figure

subplot(3,2,1)
plot(x,exacterr,'o-',x,ee,'*');
subplot(3,2,2)
plot(x,exacterr-ee,'^-')


load('err20.mat');
subplot(3,2,3)
plot(x,exacterr,'o-',x,ee,'*');
ylabel('error')
legend('exact','computed')

subplot(3,2,4)
plot(x,exacterr-ee,'^-')
ylabel('difference in error: e-$\tilde{e}$','interpreter','latex')

load('err40.mat');
subplot(3,2,5)
plot(x,exacterr,'o-',x,ee,'*');
subplot(3,2,6)
plot(x,exacterr-ee,'^-')
