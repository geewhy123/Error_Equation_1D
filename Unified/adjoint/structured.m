clear all
close all

N = 5;
xc = linspace(-0.5/N,1+0.5/N,N+2)
ue = zeros(N+2,1);
ve = zeros(N+2,1);
for i = 2:N+1
xr = xc(i)+0.5/N;
xl = xc(i)-0.5/N;
f(i) = N*pi*(cos(pi*xr)-cos(pi*xl));
ue(i) = (N/pi)*(-cos(pi*xr)+cos(pi*xl));
g(i) = (N)*(pi^3/4)*(xr^2/2-xr^3/3-xl^2/2+xl^3/3);
ve (i) = (N)*(-pi^3/48)*(xr^5/5-2*xr^4/4+xr^2/2-xl^5/5+2*xl^4/4-xl^2/2);
end
A = zeros(N+2,N+2);

for i = 2:N+1
    for j = 2:N+1
        if(i-j==1 || j-i==1)
            A(i,j) = 1;
        elseif(i==j)
            A(i,j)=-2;
        else A(i,j)=0;
        end
            
    end
end
A(2,1) = 1;
A(N+1,N+2) = 1;
A(1,1) = 1;
A(1,2) = 1;
A(N+2,N+1)=1;
A(N+2,N+2)=1;

b = [0; (1/N^2)*f(2:N+1)';0];
c = [0; (1/N^2)*g(2:N+1)';0];
u=A\b
v = A\c
plot(xc,u,'*',xc,ue)
u(1) = NaN;
u(N+2) = NaN;
ue(1) = NaN;
ue(N+2) = NaN;
v(1) = NaN;
v(N+2) = NaN;
ve(1) = NaN;
ve(N+2) = NaN;

[max(abs(ue-u)) max(abs(ve-v))]

J1 = 0;
J2 = 0;
for i = 2:N+1
   J1 = J1+ u(i)*g(i)/N;
   J2 = J2+ v(i)*f(i)/N;
end

[J1 J2]