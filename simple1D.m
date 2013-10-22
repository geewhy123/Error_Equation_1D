clear all
close all
N = 5;
h = 1/N;
A = zeros(N+2,N+2);

x = zeros(N+2,1);
for i = 1:N+2
   x(i) = (i-1.5)*h; 
end

for i = 1:N+2
   for j = 1:N+2
      if(i==j)
         A(i,i) = -2/h^2; 
      end
      if(i+1==j || i-1==j)
         A(i,j) = 1/h^2; 
      end
   end
end

A(1,1) = 1;
A(1,2) = 1;
A(N+2,N+1)=1;
A(N+2,N+2)=1;

%A(1,N) = 1;
%A(N,1) = 1;

f = zeros(N+2,1);
for i = 1:N+2
    f(i) = exp(-x(i)^2)*(-2+4*x(i)^2); 
end

f(1)=2*1;
f(N+2)=2*exp(-1);
u = A\f;



E = zeros(N+2,1);
s = zeros(N+2,1);
for i = 1:N+2
    s(i) =  %(-2+4*((i-1.5)*h)^2)*(exp(-((i-1.5)*h)^2));
    s(i) = s(i) - (1/h)*(-2*((i-1)*h)*exp(-((i-1)*h)^2)+2*((i-2)*h)*exp(-((i-2)*h)^2));
   
end
%s = -s;
s(1)=s(1)+2*(1-(u(1)+u(2))/2);
s(N+2)=s(N+2)+2*(exp(-1)-(u(N+1)+u(N+2))/2);
%s = -s;
E=A\s;


x = x(2:N+1);
u = u(2:N+1);
ue = exp(-x.^2);
plot(x,u,x,ue);

figure
plot(x,ue-u);

L2err=max(abs(ue-u))
e = ue-u;
figure 
%plot(x,u+e,x,ue);
plot(x,u,'*')
%ubar = zeroes(N,1);
%for i = 1:10*N
%    if(
%   ubar(i)  
%end


E = E(2:N+1);




figure
plot(x,E,x,e)
errorerror=max(abs(E-e))
%figure
%max(abs(s(2:N+1)-f(2:N+1)))
%plot(x,s(2:N+1),x,f(2:N+1))

%max(abs(ue-q(2:N+1)))