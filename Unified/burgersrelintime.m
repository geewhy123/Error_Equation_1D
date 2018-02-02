%bead9be63fdd706ec4beb34a895a86d9e16e4389 (HEAD)
% 02-01-18: new branch: 1a45d4454a9622433e918923186a26c446135c50
close all
clear all
set(0,'DefaultFigureVisible','off')
% N = 10:10:100;
N = [10 20 40 80];
ntimes = 5;
for i = 1:length(N);
    tic
    for j = 1:ntimes
%     [err200(i),~,~,~,~ ,~ ]  = errordriver(N(i),2,0,0,1/3,'D',0,'D',-2*tanh(1),0,'','BurgersVisc','SS',[0.2 0.2 0.2],'HC','NewtonError',0);
%     [errerr244(i),~,err(i),~,~ ,~ ]  = errordriver(N(i),2,4,4,1/3,'D',0,'D',-2*tanh(1),0,'','BurgersVisc','SS',[0.2 0.2 0.2],'HC','NewtonError',0);
% [errerr266(i),~,err(i),~,~ ,~ ]  = errordriver(N(i),2,6,6,1/3,'D',0,'D',-2*tanh(1),0,'','BurgersVisc','SS',[0.2 0.2 0.2],'HC','NewtonError',0);
%  [errerr244Relin66(i),~,err(i),~,~ ,~ ]  = errordriver(N(i),2,4,4,1/3,'D',0,'D',-2*tanh(1),0,'','BurgersVisc','SS',[0.2 0.2 0.2],'HC','NewtonError',1);

    end

    t(i) = toc;
    t(i) = t(i)/ntimes;
end
% error('1')
% r = N.^-1;
% k2 = r.^2;
% k4 = r.^4;
% k6 = r.^6;
% loglog(r,err,'-*',r,errerr244,'--o',r,errerr266,'--^',r,errerr244Relin66,'--v',r,k2,'k:',r,k4,'k-.',r,k6,'k')
% grid on
% xlabel('h')
% legend('(2,0,0)','(2,4,4) linearized','(2,6,6) linearized','(2,4,4),(2,6,6) relinearized','second order','fourth order','sixth order')
e1 = [.0264728 .005464678 .001160446 .000238275]';
e2 = [1.1379e-4 4.9037e-6 2.466e-7 1.451e-8]';
e3 = [3.921e-5 1.286e-6 8.047e-8 3.594e-9]';
e4 = [3.83e-6 3.11e-8 3.207e-10 5.867e-12]';

t1 = [0.56 0.37 0.88 2.91]';
t2 = [1.21 1.02 2.28 6.41]';
t3 = [1.30 1.16 2.45 6.71]';
t4 = [1.37 1.33 2.78 8.00]';

set(0,'DefaultFigureVisible','on')
loglog(t1,e1,t2,e2,t3,e3,t4,e4)
