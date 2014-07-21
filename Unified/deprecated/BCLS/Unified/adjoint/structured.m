clear all
close all

unif = 1/3;
N = 10;
h0 = 1/N;
x = linspace(0,1,N+1);
h = zeros(N+2,1);
for i = 2:N
    x(i) = x(i) + unif*(-1+rand*(2))*h0/3;
end
for i = 2:N+1
    h(i) = x(i)-x(i-1);
end
h(1) = h(2);
h(N+2) = h(N+1);
for i = 2:N+1
xc(i) = (x(i-1)+x(i))/2;
end
xc(1) = -xc(2);
xc(N+2) = 1+(1-xc(N+1));
ue = zeros(N+2,1);
ve = zeros(N+2,1);

for i = 2:N+1
xr = xc(i)+0.5*h(i);
xl = xc(i)-0.5*h(i);
f(i) = (1/h(i))*pi*(cos(pi*xr)-cos(pi*xl));
ue(i) = (1/h(i))*(1/pi)*(-cos(pi*xr)+cos(pi*xl));
g(i) = (1/h(i))*(pi^3/4)*(xr^2/2-xr^3/3-xl^2/2+xl^3/3);
ve (i) = (1/h(i))*(-pi^3/48)*(xr^5/5-2*xr^4/4+xr^2/2-xl^5/5+2*xl^4/4-xl^2/2);
end
A = zeros(N+2,N+2);

for i = 2:N+1
    for j = 2:N+1
%         if(i-j==1 || j-i==1)
%             A(i,j) = 1/h(i)^2;
         if(j-i==1 )
             A(i,j) = 2*(1/(h(i)*(h(i+1)+h(i))));
         elseif(i-j ==1)
             A(i,j) = 2*(1/(h(i)*(h(i-1)+h(i))));
        elseif(i==j)
%             A(i,j)=-2/h(i)^2;
            A(i,j)=-2*(1/(h(i)*(h(i+1)+h(i)))+1/(h(i)*(h(i)+h(i-1))));
        else A(i,j)=0;
        end
            
    end
end
 A(2,1) =  2*(1/(h(1)*(h(2)+h(1))));%1/h(2)^2;
A(N+1,N+2) = 1/h(N+1)^2;
A(1,1) = 1;
A(1,2) = 1;
A(N+2,N+1)=1;
A(N+2,N+2)=1;
A

b = [0; f(2:N+1)';0];
c = [0; g(2:N+1)';0];
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
   J1 = J1+ u(i)*g(i)*h(i);
   J2 = J2+ v(i)*f(i)*h(i);
end

[J1 J2]