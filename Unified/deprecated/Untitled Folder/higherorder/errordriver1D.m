close all
clear all
rng(1234);
load('20-o2')%20hru.mat')
u = lu20;
N = 20;
h0 = 1/N;
k=0.0007*(20/N)^2    ;
X = zeros(N+1,1);
for i = 1:N+1
   X(i) = (i-1)*h0; 
   if(i>1 && i < N+1)
   X(i) = X(i) + 0*(-1+rand*(2))*h0/3;
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
for i = 2:N+1
    
    f(i) = (1/h(i))*(-2*(x(i)+h(i)/2)*exp(-(x(i)+h(i)/2)^2)+2*(x(i)-h(i)/2)*exp(-(x(i)-h(i)/2)^2));
  %    f(i) = (1/h(i))*( log(1+(x(i)+h(i)/2))+((x(i)+h(i)/2)-1)/((x(i)+h(i)/2)+1) - (log(1+(x(i)-h(i)/2))+((x(i)-h(i)/2)-1)/((x(i)-h(i)/2)+1)) );%
end
global AD
AD = computepseudo(N,x,h);


ue = zeros(N+2,1);
for i = 1:N+2
ue(i) = (1/h(i))*((pi^(1/2)*erf(x(i)+h(i)/2 ))/2-(pi^(1/2)*erf(x(i)-h(i)/2 ))/2);
%ue(i) = (1/h(i))*((0.5*(1+(x(i)+h(i)/2))^2*log(1+(x(i)+h(i)/2))+7/4+1.5*(x(i)+h(i)/2)-0.25*(x(i)+h(i)/2)^2-2*(1+(x(i)+h(i)/2))*log(1+(x(i)+h(i)/2)))-((0.5*(1+(x(i)-h(i)/2))^2*log(1+(x(i)-h(i)/2))+7/4+1.5*(x(i)-h(i)/2)-0.25*(x(i)-h(i)/2)^2-2*(1+(x(i)-h(i)/2))*log(1+(x(i)-h(i)/2)))));
end


[R,uxx,Z] =computeres(u,x,h,N,f);
r=max(abs(R))





e = zeros(N+2,1);
ee =zeros(N+2,1);


tic
olderr = 1;
%esterr= zeros(,
T = 1;
for j = 1:100000

    
    %error based on u.....
[errerr(j),Z]=unstructuredrecon4(e,x,h,N,0,0);
%if(j > 50 && ((abs(errerr(j)-olderr)/olderr < 1e-15) || abs(errerr(j)-olderr)<1e-15))
if(j > 50 && k*s <1e-15)
    T = (1:1:j)*k;
    break
end
s=0;
olderr = errerr(j);
for i= 2:N+1
    y = Z(:,i);
    yr = Z(:,i+1);
    yl = Z(:,i-1);
    %need averaging flux
upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2   + 4*y(5)*(h(i)/2)^3;
upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2   + 4*yr(5)*(-h(i+1)/2)^3;
upr = (upr1+upr2)/2;
upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2   + 4*y(5)*(-h(i)/2)^3;
upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2   + 4*yl(5)*(h(i-1)/2)^3;
upl = (upl1+upl2)/2;

if i==2
    upl = upl1;
end
if i==N+1
   upr = upr1; 
end

ee(i) = e(i) + k*((upr-upl)/h(i)+R(i)); 
s=max(s,abs((upr-upl)/h(i)+R(i)));
end

e = ee;

T = (1:1:j)*k;


s
end
toc



exacterr = ue-u;
exacterr = exacterr(2:N+1);
x = x(2:N+1);
figure
plot(x,exacterr,'o-',x,ee(2:N+1),'*');
max(abs(exacterr-ee(2:N+1)))
ee = ee(2:N+1);
save('test','exacterr','ee','x')