


clear all
close all

N = [20 40 80 ];

for j = 1:length(N)
%[err(j),x,exacterr,ee] = errordriver(N(j),2,2,4,0,1000);
[perterr2(j),x,cverr2(j),ee] = errordriver(N(j),4,4,6,0,4);
[perterr3(j),x,cverr3(j),ee] = errordriver(N(j),4,4,6,0,5);
[perterr4(j),x,cverr4(j),ee] = errordriver(N(j),4,4,6,0,6);
%[perterr5(j),x,cverr5(j),ee] = errordriver(N(j),4,4,6,0,8);
end

P = N.^-4;
%loglog(N,err)
hold on
loglog(N,P,N,perterr2,N,perterr3,N,perterr4,N,perterr5);