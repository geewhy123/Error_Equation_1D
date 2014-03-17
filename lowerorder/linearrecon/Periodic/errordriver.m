%qx
function [errerr2,x,cverr2,exacterr,ee  ] = errordriver( N,p,q,r ,unif,bta,tlim,tord,physics)
%DRIVER Summary of this function goes here
%   Detailed explanation goes here

%close all

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


%u=uu;
% u(1) = NaN;
% u(N+2) = NaN;
% size(u)
% error('1')

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
    [upr,upl] = reconflux(ue,Z,f,k,h,i,N,p,physics,NaN,NaN,NaN,NaN);
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
        U(:,j) = u;
tt = k*j;
    
%[Z]=unstructuredrecon(u,x,h,N,NaN,NaN,p);


if((d*k<1e-15)||(tt>=tlim))
     
[uu,d] = updatesoln(u,x,f,k,h,N,p,tord,physics,NaN,NaN,NaN,NaN);
u = uu;
    d
    tt
    T = (1:1:j)*k;
 U(:,j+1) = u;
 nSteps = j;
 
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


[uu,d] = updatesoln(u,x,f,k,h,N,p,tord,physics,NaN,NaN,NaN,NaN);

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

figure
hold on

uder = zeros(size(U));
U(1,:) = NaN;
U(N+2,:) = NaN;
uder(1,:) = NaN;
uder(N+2,:) = NaN;
for i = 4:size(U,2)-3
w=  U(:,i);
for j = 2:N+1
%%%uder(j,i) = (U(j,i+1)-U(j,i-1))/(2*k);%(w(j+1)-w(j-1))/(2*k);
%%%%uder(j,i) = ((1/12)*U(j,i-2)+(-2/3)*U(j,i-1)+(2/3)*U(j,i+1)+(-1/12)*U(j,i+2))/(k);
uder(j,i) = ((-1/60)*U(j,i-3)+(3/20)*U(j,i-2)+(-3/4)*U(j,i-1)+(3/4)*U(j,i+1)+(-3/20)*U(j,i+2)+(1/60)*U(j,i+3))/(k);
end
end

for j = 2:N+1
%%%uder(j,1) = (-3*U(j,1)+4*U(j,2)-U(j,3))/(2*k);%(-3*w(2)+4*w(3)-w(4))/(2*k);
%%%%uder(j,1) = ((-25/12)*U(j,1)+4*U(j,2)-3*U(j,3)+(4/3)*U(j,4)-(1/4)*U(j,5))/(k);
%%%%uder(j,2) = ((-25/12)*U(j,2)+4*U(j,3)-3*U(j,4)+(4/3)*U(j,5)-(1/4)*U(j,6))/(k);
uder(j,1) = ((-49/20)*U(j,1)+6*U(j,2)-(15/2)*U(j,3)+(20/3)*U(j,4)-(15/4)*U(j,5)+(6/5)*U(j,6)-(1/6)*U(j,7))/(k);
uder(j,2) = ((-49/20)*U(j,2)+6*U(j,3)-(15/2)*U(j,4)+(20/3)*U(j,5)-(15/4)*U(j,6)+(6/5)*U(j,7)-(1/6)*U(j,8))/(k);
uder(j,3) = ((-49/20)*U(j,3)+6*U(j,4)-(15/2)*U(j,5)+(20/3)*U(j,6)-(15/4)*U(j,7)+(6/5)*U(j,8)-(1/6)*U(j,9))/(k);

%%%uder(j,size(U,2)) = (-3*U(j,size(U,2))+4*U(j,size(U,2)-1)-U(j,size(U,2)-2))/(2*k);%(-3*w(N+1)+4*w(N)-w(N-1))/(2*k);
%%%%uder(j,size(U,2)-1) = -((-25/12)*U(j,size(U,2)-1)+4*U(j,size(U,2)-2)-3*U(j,size(U,2)-3)+(4/3)*U(j,size(U,2)-4)-(1/4)*U(j,size(U,2)-5))/(k);
%%%%uder(j,size(U,2)) = -((-25/12)*U(j,size(U,2))+4*U(j,size(U,2)-1)-3*U(j,size(U,2)-2)+(4/3)*U(j,size(U,2)-3)-(1/4)*U(j,size(U,2)-4))/(k);
uder(j,size(U,2)-2) = -((-49/20)*U(j,size(U,2)-2)+6*U(j,size(U,2)-3)-(15/2)*U(j,size(U,2)-4)+(20/3)*U(j,size(U,2)-5)-(15/4)*U(j,size(U,2)-6)+(6/5)*U(j,size(U,2)-7)-(1/6)*U(j,size(U,2)-8))/(k);
uder(j,size(U,2)-1) = -((-49/20)*U(j,size(U,2)-1)+6*U(j,size(U,2)-2)-(15/2)*U(j,size(U,2)-3)+(20/3)*U(j,size(U,2)-4)-(15/4)*U(j,size(U,2)-5)+(6/5)*U(j,size(U,2)-6)-(1/6)*U(j,size(U,2)-7))/(k);
uder(j,size(U,2))   = -((-49/20)*U(j,size(U,2))+6*U(j,size(U,2)-1)-(15/2)*U(j,size(U,2)-2)+(20/3)*U(j,size(U,2)-3)-(15/4)*U(j,size(U,2)-4)+(6/5)*U(j,size(U,2)-5)-(1/6)*U(j,size(U,2)-6))/(k);
end
uder(:,:);
size(U)
size(uder)

%plot(x,(U(:,end)-U(:,end-2))/(2*k),'o',x,uder(:,end-1),'*')
plot(x,uder(:,end),'*')
hold on

%XX = linspace(0,1,100);
%YY = 2*pi*cos(2*pi*(XX+1));
%YY = -400*pi^2*exp(-0.4*pi^2)*sin(2*pi*XX);
XX = x;
YY = zeros(size(XX));
for i = 2:N+1
YY(i) = (1/h(i))*(sin(2*pi*(XX(i)+h(i)/2+1))-sin(2*pi*(XX(i)-h(i)/2+1)));
end
plot(XX,YY,'+')
max(abs(YY-uder(:,end)))

k
uder(:,end)

T=(0:1:nSteps)*k;
%cs =csapi(T,U(3,:))


for j = 2:N+1
sp = spapi(6,T,U(j,:));
gsp(j) = sp;
%figure
fnplt(sp)
%hold on
%plot(T,U(3,:),'*')
end
% 
% M = 40;
% xx = linspace(0,1,M);
% yy = exp(sin(pi*xx));
% sp = spapi(6,xx,yy);
% figure 
% fnplt(sp)
% fnval(fnder(sp),.25)-exp(sin(pi*.25))*pi*cos(pi*.25)
% fnval(sp,.5)-exp(1)
%  error('1')
% size(sp.coefs)
% gsp
% figure
% fnplt(sp)
% fnval(sp,.3)
% error('1')

% % % [m,n] = size(uder);
% % % UT = zeros(m,n);
% % % for i = 2:N+1
% % %     sp = gsp(i);
% % %     for j = 1:n
% % %         ttt = (j-1)*k;
% % %       UT(i,j) = fnval(fnder(sp),ttt);     
% % %     end
% % % end
% % % norm(UT(2:N+1,:)-uder(2:N+1,:));
% % % %UT
% % % %uder



if(q>0 && r > 0)
    
    clearvars -except u N p q r unif FI bta f cverr2 v k ue u0 tlim tord uo physics uder nSteps gsp U
    
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

global TEND
TEND = tlim;

global UU
UU = U;


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



% [R0,uxx,Z] =computeres(u,x,k,h,N,f,r,physics,uder,1);
% res=max(abs(R0))
% R0;
uold = u;
u=u0;

global R
tt=0;
  [R(:,1),uxx,Z] =computeres(u,x,k,h,N,f,r,physics,uder,1,0,gsp);
for j = 1:nSteps
%     uder(:,j) = (1/h(i))*(sin(2*pi*(x(i)+h(i)/2+(j-1)*k))-sin(2*pi*(x(i)-h(i)/2+(j-1)*k)));
%     for i = 2:N+1
%     u(i) =  (1/h(i))*(-1/(2*pi))*(cos(2*pi*(x(i)+h(i)/2+(j-1)*k)) -cos(2*pi*(x(i)-h(i)/2+(j-1)*k)));
%     end

AD = computepseudo(N,x,h,p);    
     [u,d] = updatesoln(u,x,f,k,h,N,p,tord,physics,NaN,NaN,NaN,NaN);%uder,j,tt,gsp);
    %%recon(p)??every time step

    
 AD = computepseudo(N,x,h,r);
   [R(:,j+1),uxx,Z] =computeres(U(:,j+1),x,k,h,N,f,r,physics,uder,j+1,tt+k,gsp);
   
tt = tt+k;
   
%    if(abs(j*k-0.2) < 1e-5)
%       max(abs(R(:,j)))
%       error('1')
%    end
end


%R(:,1)=0;
%uder
%R
max(abs(R(:,end)))

%error('1')
%uder(:,end)
% error('1')

% R0-R(:,end)
% error('1')


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
% % % uu = ue;
% % % uu=uu+h0^bta*v;
% % % AD = computepseudo(N,x,h,q);
% % % FIq = zeros(N+2,1);
% % %     [Z]=unstructuredrecon(uu,x,h,N,NaN,NaN,q);
% % % for i = 2:N+1
% % %     [upr,upl] = reconflux(uu,Z,f,k,h,i,N,q,physics);
% % %     FIq(i) = (upr-upl)/h(i)-f(i);
% % % end
% % % FIq(N+2) = NaN;  
% % % 
% % % AD = computepseudo(N,x,h,p);
% % % FIp = zeros(N+2,1);
% % %     [Z]=unstructuredrecon(uu,x,h,N,NaN,NaN,p);
% % % for i = 2:N+1
% % %     [upr,upl] = reconflux(uu,Z,f,k,h,i,N,p,physics);
% % %     FIp(i) = (upr-upl)/h(i)-f(i);
% % % end
% % % FIp(N+2) = NaN;  
% % % 
% % % 
% % % 
% % % FI-FIq;

%R = -FInew;
%R = -FIp;
%R=(-FIp-FIq)/2
% v = rand(N+2,1);
% v = v./norm(v);
% R = R+h0^bta*v'
% SS = dot(h(2:N+1),R(2:N+1))
% R=R-SS
%%%%%
% 
% for j = 2:N+1
% sp = spapi(6,T,R(j,:));
% Rsp(j) = sp;
% 
% end


global M
M = nSteps;


 AD = computepseudo(N,x,h,q);

 
 
T = 1;
s=1;

for j = 1:100000


%  AD = computepseudo(N,x,h,r);
% [R,uxx,Z] =computeres(u,x,k,h,N,f,r,physics,uder,j);
%  AD = computepseudo(N,x,h,q);
E(:,j) = e;

    TT = k*j;
    
%%%%[Z]=unstructuredrecon(e,x,h,N,NaN,NaN,q);



if( ((s*k<1e-15)||(TT>=tlim)) || (j > nSteps))
% AD = computepseudo(N,x,h,r);
% [R,uxx,Z] =computeres(u,x,k,h,N,f,r,physics,uder,j);
%  AD = computepseudo(N,x,h,q);
[ee,s] = updatesoln(e,x,-R(:,j),k,h,N,q,tord,physics,NaN,NaN,TT,gsp);



e = ee;
    s
    TT
    T = (1:1:j)*k;
j
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

[ee,s] = updatesoln(e,x,-R(:,j),k,h,N,q,tord,physics,NaN,NaN,TT,gsp);


e = ee;

T = (1:1:j)*k;


if(mod(j,100)==0)
s
end
end




exacterr = ue-u;

norm(u(2:N+1))

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




clear global
end
