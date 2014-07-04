clear all
close all
rng(1234);
p=6
q = 1+1e-9
N = 20
x = [rand(1,N) 0 1];

y = exp(-x.^2);

sp = spapi(p,x,y);


abs(fnval(sp,q)-exp(-q^2))
plot(x,y,'o')