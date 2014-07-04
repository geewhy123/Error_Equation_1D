close all

load('N5.mat')
load('N10.mat')
load('N20.mat')
load('N40.mat')
T40 = T;
err40=err;
u40 = u;
load('N80.mat')


semilogy(T5,err5,'*',T10,err10,'o',T20,err20,'s',T40,err40,'^',T,err,'v')

err5(end)/err10(end)
err10(end)/err20(end)
err20(end)/err40(end)
err40(end)/err(end)

min(err5)/min(err10)
min(err10)/min(err20)
min(err20)/min(err40)
min(err40)/min(err)

set(gca,'FontSize',18)
xlabel('t')
h=ylabel('error,sup $$|u_i(x)$$-$$\tilde{u}_i(x)|$$')
legend('N=5','N=10','N=20','N=40','N=80');
set(h,'Interpreter','latex','fontsize',18)

figure
clear all
load('N5u.mat')
load('N10u.mat')
load('N20u.mat')
semilogy(T5,err5,'*',T10u,err10u,'o',T20u,err20u,'s')
ylim([1e-8 1e2])

err5(end)/err10u(end)
err10u(end)/err20u(end)
