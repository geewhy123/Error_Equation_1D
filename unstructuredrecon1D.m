clear all
close all
N = 10;
rng(1234);
h0 = 1/N;
X = zeros(N+1,1);
for i = 1:N+1
   X(i) = (i-1)*h0; 
   if(i>1 && i < N+1)
   X(i) = X(i) + (-1+rand*(2))*h0/3;
   end
end
x = zeros(N+2,1);
for i = 2:N+1
    x(i) = (X(i-1)+X(i))/2;
end
x(1) = -x(2);
x(N+2) = 1+(1-x(N+1));
%plot(X,1,'*',x,1,'o')
%xlim([0 1])
h = zeros(N+2,1);
for i = 2:N+1
   h(i) = X(i)-X(i-1); 
end

h(1) = h(2);
h(N+2) = h(N+1);



f = zeros(N+2,1);
for i = 1:N+2
    if i==1
       f(i) = 0; 
    else
    f(i) = (1/h(i))*(-2*(x(i)+h(i)/2)*exp(-(x(i)+h(i)/2)^2)+2*(x(i-1)+h(i-1)/2)*exp(-(x(i-1)+h(i-1)/2)^2));
    end
end

f(1)=2*1;
f(N+2)=2*exp(-1);

plot(x,f,'o')