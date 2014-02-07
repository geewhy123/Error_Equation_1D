%q
function [errerr2,x,exacterr,ee  ] = errordriver( N,p,q,r ,unif)
%DRIVER Summary of this function goes here
%   Detailed explanation goes here

close all

if(p>0)
%N = 80;
rng(1234);
h0 = 1/N;
k = .004*(10/N)^2;%0.00006;
X = zeros(N+1,1);
for i = 1:N+1
   X(i) = (i-1)*h0; 
   if(i>1 && i < N+1)
   X(i) = X(i) + unif*(-1+rand*(2))*h0/3;
   end
end
x = zeros(N+2,1);
for i = 2:N+1
    x(i) = (X(i-1)+X(i))/2;
end

x(1) = 0-(1-x(N+1));%-x(2);
x(N+2) = 1+x(2);%1+(1-x(N+1));

h = zeros(N+2,1);
for i = 2:N+1
   h(i) = X(i)-X(i-1); 
end

h(1) = h(N+1);
h(N+2) = h(2);

u = zeros(N+2,1);

f = zeros(N+2,1);
for i = 1:N+2
    if i==1
   ;   %%% f(i) = 0; 
    else
    %f(i) = (1/h(i))*(-2*(x(i)+h(i)/2)*exp(-(x(i)+h(i)/2)^2)+2*(x(i-1)+h(i-1)/2)*exp(-(x(i-1)+h(i-1)/2)^2));
    xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;
    %f(i)= (1/h(i))*(1/(12*pi))*(3*cos(2*pi*xr)-cos(6*pi*xr)-3*cos(2*pi*xl)+cos(6*pi*xl));
    %f(i) = (1/h(i))*(-2*pi);%*(exp(cos(2*pi*xr))*   (2*(cos(pi*xr)^2-    4*(cos(pi*xr))^4+1)));% - (exp(cos(2*pi*xl))*(2*(cos(pi*xl))^2-4*(cos(pi*xl))^4+1));
    %f(i)=(1/h(i))*((-1/pi)*sin(2*pi*xr)+(1/(pi^3))*(sin(4*pi*xr))+(1/pi)*(sin(2*pi*xl))-(1/(pi^3))*(sin(4*pi*xl)));
   
f(i) = (1/h(i))*(-2*pi)*(sin(2*pi*xr)-sin(2*pi*xl));
%f(i) =     (1/h(i))*((-4*pi^2-)/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl));


   % f(i) = (1/h(i))*(-4*pi^2)*( (exp(1)^3*sin(2*pi*xr)+1)/(sin(2*pi*xr)+exp(1)^3)^2 - (exp(1)^3*sin(2*pi*xl)+1)/(sin(2*pi*xl)+exp(1)^3)^2);
   
    end
    f(1) = f(N+1);
    %u(i) = (1/h(i))*((x(i)+h(i)/2)^4-((x(i)-h(i)/2)^4))/4;%exp(-(x(i)-0.5)^2);
    %u(i) = (1/h(i))*((x(i)+h(i)/2)^5-((x(i)-h(i)/2)^5))/5;%exp(-(x(i)-0.5)^2);
    %u(i) = (1/h(i))*((pi^(1/2)*erf(x(i)+h(i)/2 ))/2-(pi^(1/2)*erf(x(i)-h(i)/2 ))/2);
    %u(i) = (1/h(i))*(-1/(2*pi))*(exp(cos(2*pi*(x(i)+h(i)/2)))-exp(cos(2*pi*(x(i)-h(i)/2))));    
    %u(i) = (1/h(i))*(-pi)*(cos(2*pi*(x(i)+h(i)/2))-3*cos(6*pi*(x(i)+h(i)/2))+cos(2*pi*(x(i)-h(i)/2))+3*cos(6*pi*(x(i)-h(i)/2)));
      xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;
   
    %u(i) = (1/h(i))*(4*pi*sin(2*pi*xr)-(16/pi)*sin(4*pi*xr)-4*pi*sin(2*pi*xl)+(16/pi)*sin(4*pi*xl));
    
   u(i) = (1/h(i))*(1/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl));
  %u(i) = (1/h(i))*(log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl)));
end

 u(1) = NaN;
 u(N+2) = NaN;
ue = u;
 u
 f
%u(:)=0;

%f(1)=2*1;
%f(N+2)=2*exp(-1);


%hold on

uu = zeros(N+2,1);

global AD
AD = computepseudo(N,x,h,p);
%AD
%cond(AD(:,:,4))
%%%plot(x,u,'*')

    [err,Z]=unstructuredrecon(u,x,h,N,NaN,NaN,p);
for i = 2:N+1
    [upr,upl] = reconflux(u,Z,f,k,h,i,N,p);
    FI(i) = (upr-upl)/h(i)-f(i);
    d2u(i) = (upr-upl)/h(i);
end
d2u
FI
fi= max(abs(FI))
%Z
%Z
%error('1')


d=1;
tic
T = 1;
for j = 1:100000

    
[err(j),Z]=unstructuredrecon(u,x,h,N,NaN,NaN,p);



%if(j > 10 && ((abs(err(j)-olderr)/olderr < 1e-8) || abs(err(j)-olderr)<1e-15))
if(d*k<1e-15)
    T = (1:1:j)*k;
    break
end

d=0;



for i= 2:N+1
    
    [upr,upl] = reconflux(u,Z,f,k,h,i,N,p);
    
uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
d = max(d,abs((upr-upl)/h(i)-f(i)));




end

u = uu;

T = (1:1:j)*k;

if(mod(j,100)==0)
    d
end

end
toc

save('test','u','T','err','x')


cverr1 = sum(abs(ue(2:N+1)-u(2:N+1)))/N
cverr2 = sqrt(sum((ue(2:N+1)-u(2:N+1)).^2)/N)

cverrinf=max(abs(ue-u))


%m = N/2;
%FI = (-u(m-2)+16*u(m-1)-30*u(m)+16*u(m+1)-u(m+2)) /(12*h(m)^2)-f(m)

u(1) = NaN;
u(N+2) = NaN;
%%%plot(x,u,'*',x,ue,'o')
%%%figure
%%%plot(x,ue-u,'x')

end

%h
%x




if(q>0 && r > 0)
    
    clearvars -except u N p q r unif
    
%Error equation
%clear all
rng(1234);
load('test.mat')

h0 = 1/N;
k=0.0008  *(20/N)^2  ;
X = zeros(N+1,1);
for i = 1:N+1
   X(i) = (i-1)*h0; 
   if(i>1 && i < N+1)
   X(i) = X(i) + unif*(-1+rand*(2))*h0/3;
   end
end
x = zeros(N+2,1);
for i = 2:N+1
    x(i) = (X(i-1)+X(i))/2;
end
x(1) = 0-(1-x(N+1));%-x(2);
x(N+2) = 1+x(2);%1+(1-x(N+1));
%plot(X,1,'*',x,1,'o')
%xlim([0 1])
h = zeros(N+2,1);
for i = 2:N+1
   h(i) = X(i)-X(i-1); 
end

h(1) = h(N+1);
h(N+2) = h(2);


f = zeros(N+2,1);
for i = 2:N+1
   
    xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;
    %f(i)= (1/h(i))*(1/(12*pi))*(3*cos(2*pi*xr)-cos(6*pi*xr)-3*cos(2*pi*xl)+cos(6*pi*xl));
    %f(i) = (1/h(i))*(-2*pi);%*(exp(cos(2*pi*xr))*   (2*(cos(pi*xr)^2-    4*(cos(pi*xr))^4+1)));% - (exp(cos(2*pi*xl))*(2*(cos(pi*xl))^2-4*(cos(pi*xl))^4+1));
    %f(i)=(1/h(i))*((-1/pi)*sin(2*pi*xr)+(1/(pi^3))*(sin(4*pi*xr))+(1/pi)*(sin(2*pi*xl))-(1/(pi^3))*(sin(4*pi*xl)));
  
%f(i) = (1/h(i))*(-4*pi^2)*( (exp(1)^3*sin(2*pi*xr)+1)/(sin(2*pi*xr)+exp(1)^3)^2 - (exp(1)^3*sin(2*pi*xl)+1)/(sin(2*pi*xl)+exp(1)^3)^2); 
f(i) = (1/h(i))*(-2*pi)*(sin(2*pi*xr)-sin(2*pi*xl));


% f(i) = (1/h(i))*(-2*(x(i)+h(i)/2)*exp(-(x(i)+h(i)/2)^2)+2*(x(i-1)+h(i-1)/2)*exp(-(x(i-1)+h(i-1)/2)^2));
  %    f(i) = (1/h(i))*( log(1+(x(i)+h(i)/2))+((x(i)+h(i)/2)-1)/((x(i)+h(i)/2)+1) - (log(1+(x(i)-h(i)/2))+((x(i)-h(i)/2)-1)/((x(i)-h(i)/2)+1)) );%
end
% global PS3
% PS3 = computepseudo3(N,x,h);
% global PS1
% PS1 = computepseudo(N,x,h);
% global PS4
% PS4 = computepseudo4(N,x,h);
 global AD
 AD = computepseudo(N,x,h,r);

ue = zeros(N+2,1);
for i = 1:N+2
        xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;
   
    %u(i) = (1/h(i))*(4*pi*sin(2*pi*xr)-(16/pi)*sin(4*pi*xr)-4*pi*sin(2*pi*xl)+(16/pi)*sin(4*pi*xl));
  
    %ue(i) = (1/h(i))*(log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl)));
  ue(i) = (1/h(i))*(1/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl));
  
%ue(i) = (1/h(i))*((pi^(1/2)*erf(x(i)+h(i)/2 ))/2-(pi^(1/2)*erf(x(i)-h(i)/2 ))/2);
%ue(i) = (1/h(i))*((0.5*(1+(x(i)+h(i)/2))^2*log(1+(x(i)+h(i)/2))+7/4+1.5*(x(i)+h(i)/2)-0.25*(x(i)+h(i)/2)^2-2*(1+(x(i)+h(i)/2))*log(1+(x(i)+h(i)/2)))-((0.5*(1+(x(i)-h(i)/2))^2*log(1+(x(i)-h(i)/2))+7/4+1.5*(x(i)-h(i)/2)-0.25*(x(i)-h(i)/2)^2-2*(1+(x(i)-h(i)/2))*log(1+(x(i)-h(i)/2)))));
end


[R,uxx,Z] =computeres(u,x,h,N,f,r);
res=max(abs(R))

%error('1')


e = zeros(N+2,1);
ee =zeros(N+2,1);




 AD = computepseudo(N,x,h,q);

tic
olderr = 1;
%esterr= zeros(,
T = 1;
for j = 1:100000

    
    
    %error based on u.....
[errerr(j),Z]=unstructuredrecon(e,x,h,N,NaN,NaN,q);
%if(j > 50 && ((abs(errerr(j)-olderr)/olderr < 1e-15) || abs(errerr(j)-olderr)<1e-15))
if(j > 50 && k*s <1e-15)
    T = (1:1:j)*k;
    break
end
s=0;
olderr = errerr(j);
for i= 2:N+1


    
    [upr,upl] = reconflux(e,Z,-R,k,h,i,N,q);



ee(i) = e(i) + k*((upr-upl)/h(i)+R(i)); 
s=max(s,abs((upr-upl)/h(i)+R(i)));
end

e = ee;

T = (1:1:j)*k;


if(mod(j,100)==0)
s
end
end
toc



exacterr = ue-u;
exacterr = exacterr(2:N+1);
x = x(2:N+1);
figure
plot(x,exacterr,'o-',x,ee(2:N+1),'*');
errerr= max(abs(exacterr-ee(2:N+1)));



ee = ee(2:N+1);
ue=ue(2:N+1);


errerr1 = sum(abs(exacterr-ee))/N
errerr2 = sqrt(sum((exacterr-ee).^2)/N)

errerrinf=max(abs(exacterr-ee))


save('t','exacterr','ee','x')



end

end
