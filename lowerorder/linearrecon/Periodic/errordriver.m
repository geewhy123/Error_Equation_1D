%qx
function [errerr2,x,cverr2,exacterr,ee  ] = errordriver( N,p,q,r ,unif,bta,tlim,tord,physics)
%DRIVER Summary of this function goes here
%   Detailed explanation goes here

close all

if(p>0)

rng(1234);
h0 = 1/N;

CFL = 0.4;
k = CFL*h0;

if(strcmp(physics,'Poisson')==1)
    k = k*h0;
end

%k = .004*(10/N)^2
%k=k/4
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

% ue = zeros(N+2,1);
% u0 = zeros(N+2,1);
% f = zeros(N+2,1);
% for i = 2:N+1
%     
%     xl = x(i)-h(i)/2;
%     xr = x(i)+h(i)/2;
%     %f(i)= (1/h(i))*(1/(12*pi))*(3*cos(2*pi*xr)-cos(6*pi*xr)-3*cos(2*pi*xl)+cos(6*pi*xl));
%     %f(i) = (1/h(i))*(-2*pi);%*(exp(cos(2*pi*xr))*   (2*(cos(pi*xr)^2-    4*(cos(pi*xr))^4+1)));% - (exp(cos(2*pi*xl))*(2*(cos(pi*xl))^2-4*(cos(pi*xl))^4+1));
%     %f(i)=(1/h(i))*((-1/pi)*sin(2*pi*xr)+(1/(pi^3))*(sin(4*pi*xr))+(1/pi)*(sin(2*pi*xl))-(1/(pi^3))*(sin(4*pi*xl)));
%    
%     %%%%%f(i) = (1/h(i))*(-2*pi)*(sin(2*pi*xr)-sin(2*pi*xl));
%     %f(i) =     (1/h(i))*((-4*pi^2-)/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl));
% 
% 
%      %this%%%   f(i) = (1/h(i))*(-4*pi^2)*( (exp(1)^3*sin(2*pi*xr)+1)/(sin(2*pi*xr)+exp(1)^3)^2 - (exp(1)^3*sin(2*pi*xl)+1)/(sin(2*pi*xl)+exp(1)^3)^2);
%    %f(i) = -(1/h(i))*(2*pi)*(sin(2*pi*xr)-sin(2*pi*xl));
%     
%    f(i) = 0;
%     
%     %u(i) = (1/h(i))*((x(i)+h(i)/2)^4-((x(i)-h(i)/2)^4))/4;%exp(-(x(i)-0.5)^2);
%     %u(i) = (1/h(i))*((x(i)+h(i)/2)^5-((x(i)-h(i)/2)^5))/5;%exp(-(x(i)-0.5)^2);
%     %u(i) = (1/h(i))*((pi^(1/2)*erf(x(i)+h(i)/2 ))/2-(pi^(1/2)*erf(x(i)-h(i)/2 ))/2);
%     %u(i) = (1/h(i))*(-1/(2*pi))*(exp(cos(2*pi*(x(i)+h(i)/2)))-exp(cos(2*pi*(x(i)-h(i)/2))));    
%     %u(i) = (1/h(i))*(-pi)*(cos(2*pi*(x(i)+h(i)/2))-3*cos(6*pi*(x(i)+h(i)/2))+cos(2*pi*(x(i)-h(i)/2))+3*cos(6*pi*(x(i)-h(i)/2)));
%     xl = x(i)-h(i)/2;
%     xr = x(i)+h(i)/2;
%    
%     %u(i) = (1/h(i))*(4*pi*sin(2*pi*xr)-(16/pi)*sin(4*pi*xr)-4*pi*sin(2*pi*xl)+(16/pi)*sin(4*pi*xl));
%     
% %%%%%    ue(i) = (1/h(i))*(1/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl));
% 
% 
%  %this%% ue(i) = (1/h(i))*(log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl)));
%  
% %%%this ue(i)= (1/h(i))*((-1/(2*pi))*(100*exp(-4*pi^2*tlim))*(cos(2*pi*xr)-cos(2*pi*xl))+  (log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl))));%(1/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl)));
% ue(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*(xr+tlim)) -cos(2*pi*(xl+tlim)));
% 
% 
% %%%burgers
% 
% 
% 
%  %initial
% %this u(i) = (1/h(i))*((-1/(2*pi))*(100*(cos(2*pi*xr)-cos(2*pi*xl))) +(log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl))));%(1/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl)));
%  
% u0(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
% 
% end
% f(1) = NaN;
% f(N+2) = NaN;
%  ue(1) = NaN;
%  ue(N+2) = NaN;
%  u0(1) = NaN;
%  u0(N+2)= NaN;
 

 [u0,ue,f]=initializeexact(physics,N,x,h,tlim);
 
%u = ue;
u=u0;
uu = zeros(N+2,1);




v = rand(N+2,1);
v = v./norm(v);
%R = R+h0^bta*v'
%SS = dot(h(2:N+1),R(2:N+1))
%R=R-SS

% plot(x,u0-ue)
% max(abs(u0-ue))
% error('1')

u=u+h0^bta*v;


global AD
AD = computepseudo(N,x,h,p);
FI = zeros(N+2,1);
d2u = zeros(N+2,1);

    [Z]=unstructuredrecon(ue,x,h,N,NaN,NaN,p);
   % [Ze]= unstructuredrecon(ue,x,h,N,NaN,NaN,p);
for i = 2:N+1
    [upr,upl] = reconflux(ue,Z,f,k,h,i,N,p);
  %  [upre,uple]=reconflux(ue,Ze,f,k,h,i,N,p);
    FI(i) = (upr-upl)/h(i)-f(i);
 %   FIe(i)=(upre-uple)/h(i)-f(i);
    d2u(i) = (upr-upl)/h(i);
end


d2u;
FI(N+2) = NaN;  
FI;

%norm(FIe(2:N+1)-FI(2:N+1))

%error('1')


fi= max(abs(FI))
%Z

d=1;

T = 1;
for j = 1:100000
        
tt = k*j;
    
%[Z]=unstructuredrecon(u,x,h,N,NaN,NaN,p);


if((d*k<1e-15)||(tt>=tlim))
    
[uu,d] = updatesoln(u,x,f,k,h,N,p,tord);
u = uu;
    d
    tt
    T = (1:1:j)*k;
    break
end

d=0;

%for i= 2:N+1
    
   % [upr,upl] = reconflux(u,Z,f,k,h,i,N,p);

    
%    [uu(i)] = updatesol(
    
%uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 


%d = max(d,abs((upr-upl)/h(i)-f(i)));


%[delt]= updatesol(u,Z,f,k,h,i,N,p);
%uu(i) = u(i)+k*delt;
%d = max(d,abs(delt));


[uu,d] = updatesoln(u,x,f,k,h,N,p,tord);

%end
uo = u;
u = uu;

T = (1:1:j)*k;

if(mod(j,100)==0)
    d
end

end


save('test','u','T','x')


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

plot(x,u,x,ue)

end

T(end)

if(q>0 && r > 0)
    
    clearvars -except u N p q r unif FI bta f cverr2 v k ue u0 tlim tord uo
    
%Error equation
%clear all
rng(1234);
load('test.mat')

h0 = 1/N;
%k=0.0008  *(20/N)^2  ;
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



 global AD
 AD = computepseudo(N,x,h,r);

%ue = zeros(N+2,1);
for i = 1:N+2
        xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;
   
    %u(i) = (1/h(i))*(4*pi*sin(2*pi*xr)-(16/pi)*sin(4*pi*xr)-4*pi*sin(2*pi*xl)+(16/pi)*sin(4*pi*xl));
  
   % ue(i) = (1/h(i))*(log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl)));
  %%%%%ue(i) = (1/h(i))*(1/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl));
  
%ue(i) = (1/h(i))*((pi^(1/2)*erf(x(i)+h(i)/2 ))/2-(pi^(1/2)*erf(x(i)-h(i)/2 ))/2);
%ue(i) = (1/h(i))*((0.5*(1+(x(i)+h(i)/2))^2*log(1+(x(i)+h(i)/2))+7/4+1.5*(x(i)+h(i)/2)-0.25*(x(i)+h(i)/2)^2-2*(1+(x(i)+h(i)/2))*log(1+(x(i)+h(i)/2)))-((0.5*(1+(x(i)-h(i)/2))^2*log(1+(x(i)-h(i)/2))+7/4+1.5*(x(i)-h(i)/2)-0.25*(x(i)-h(i)/2)^2-2*(1+(x(i)-h(i)/2))*log(1+(x(i)-h(i)/2)))));
end



[R,uxx,Z] =computeres(u,x,k,h,N,f,r);
res=max(abs(R))
R
%error('1')

% error('1')


e = zeros(N+2,1);
ee =zeros(N+2,1);


% AD = computepseudo(N,x,h,q);
% FInew = zeros(N+2,1);
%     [Z]=unstructuredrecon(ue,x,h,N,NaN,NaN,q);
% for i = 2:N+1
%     [upr,upl,FInew(i)] = reconflux(ue,Z,f,k,h,i,N,q);
%     %FIq(i) = (upr-upl)/h(i)-f(i);
% end
% FIq(N+2) = NaN; 




%%%%perturb
uu = ue;
uu=uu+h0^bta*v;
AD = computepseudo(N,x,h,q);
FIq = zeros(N+2,1);
    [Z]=unstructuredrecon(uu,x,h,N,NaN,NaN,q);
for i = 2:N+1
    [upr,upl] = reconflux(uu,Z,f,k,h,i,N,q);
    FIq(i) = (upr-upl)/h(i)-f(i);
end
FIq(N+2) = NaN;  

AD = computepseudo(N,x,h,p);
FIp = zeros(N+2,1);
    [Z]=unstructuredrecon(uu,x,h,N,NaN,NaN,p);
for i = 2:N+1
    [upr,upl] = reconflux(uu,Z,f,k,h,i,N,p);
    FIp(i) = (upr-upl)/h(i)-f(i);
end
FIp(N+2) = NaN;  



FI-FIq;

%R = -FInew;
%R = -FIp;
%R=(-FIp-FIq)/2
% v = rand(N+2,1);
% v = v./norm(v);
% R = R+h0^bta*v'
% SS = dot(h(2:N+1),R(2:N+1))
% R=R-SS
%%%%%


 AD = computepseudo(N,x,h,q);


T = 1;
s=1;
for j = 1:100000
    TT = k*j;
    
[Z]=unstructuredrecon(e,x,h,N,NaN,NaN,q);



if( ((s*k<1e-15)||(TT>=tlim)))
    
[ee,s] = updatesoln(e,x,-R,k,h,N,q,tord);

e = ee;
    s
    TT
   
    T = (1:1:j)*k;
    break
end

% if(j > 50 && k*s <1e-15)
%     T = (1:1:j)*k;
%     break
% end
s=0;
%%%for i= 2:N+1
   
  %%%  [upr,upl] = reconflux(e,Z,-R,k,h,i,N,q);
% 
%     
% ee(i) = e(i) + k*((upr-upl)/h(i)+R(i)); 
% s=max(s,abs((upr-upl)/h(i)+R(i)));
% 
% 
% end

[ee,s] = updatesoln(e,x,-R,k,h,N,q,tord);


e = ee;

T = (1:1:j)*k;


if(mod(j,100)==0)
s
end
end




exacterr = ue-u;
exacterr = exacterr(2:N+1);
x = x(2:N+1);
figure
plot(x,exacterr,'o-',x,ee(2:N+1),'*');


ee = ee(2:N+1);
ue=ue(2:N+1);


errerr1 = sum(abs(exacterr-ee))/N
errerr2 = sqrt(sum((exacterr-ee).^2)/N)

errerrinf=max(abs(exacterr-ee))


figure
plot(x,exacterr-ee,'*-')

%norm(exacterr)
%[exacterr ee (exacterr-ee)*1e6]


save('t','exacterr','ee','x')


end

end
