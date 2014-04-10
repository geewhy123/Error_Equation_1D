
function [errerr2,x,cverr2,exacterr,ee  ] = errordriver( N,p,q,r ,unif,bta,tlim,tord,physics,goal)
%DRIVER Summary of this function goes here
%   Detailed explanation goes here

close all

if(p>0)

%rng(1234);
w = 0;
h0 = 1/N;

CFL = 0.4;
k = CFL*h0;

if(strcmp(physics,'Poisson')==1)
    k = k*h0;
%     if (N>10)
%     k=4*k;
%     elseif (N>20)
%     k = 8*k;
%     end
% if(strcmp(goal,'SS')==1)
    k = 2*k;
if( tlim > 2)
    k = 1.5*k;
end
end


X = zeros(N+1,1);
for i = 1:N+1
   X(i) = (i-1)*h0; 
   if(i>1 && i < N+1)
   X(i) = X(i) + unif*(-1+rand*(2))*h0/3;%0.001*sin(2*pi*X(i));%
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


 [u0,ue,f]=initializeexact(physics,N,x,h,tlim);
 

%u = ue;
u=u0;
%uu = zeros(N+2,1);



global AD
AD = computepseudo(N,x,h,p);    


d=1;

T = 1;
for j = 1:100000
        U(:,j) = u;
tt = k*j;
    


if((max(d)*k<1e-15)||(tt>=tlim))
     
[uu,d] = update('solution',u,x,f,k,h,N,p,tord,physics,NaN,NaN);
u = uu;
    max(d)
    tt
    T = (1:1:j)*k;
 U(:,j+1) = u;
 nSteps = j;
 
    break
end

d=0;


[uu,d] = update('solution',u,x,f,k,h,N,p,tord,physics,NaN,NaN);


u = uu;

T = (1:1:j)*k;

if(mod(j,100)==0)
    max(d)
end

end



cverr1 = sum(abs(ue(2:N+1)-u(2:N+1)))/N
cverr2 = sqrt(sum((ue(2:N+1)-u(2:N+1)).^2)/N)
cverrinf=max(abs(ue-u))



u(1) = NaN;
u(N+2) = NaN;
plot(x,u,'*',x,ue,'o')
figure
plot(x,ue-u,'x')


end

T(end)
assert((abs(T(end)-tlim)/tlim < 1e-4) || (strcmp(physics,'Poisson')==1 && tlim/T(end) > 2 ) ) 


global dUdt
dUdt = diffU(U,k);
gsp = NaN;

% if not using FD, then use this to spline
% % % T=(0:1:nSteps)*k;
% % % for j = 2:N+1
% % % sp = spapi(6,T,U(j,:));
% % % gsp(j) = sp;
% % % %figure
% % % %fnplt(sp)
% % % %hold on
% % % %plot(T,U(3,:),'*')
% % % end




nSteps
% 
% if(nSteps < 20)
%    error('1') 
% end

if(q>0 && r > 0)
    
    clearvars -except u N p q r unif FI bta f cverr2 v k ue u0 tlim tord uo physics uder nSteps gsp U h x goal
    
    
    
figure
hold on

global KK
KK = k;
global xx
xx = x;
global TEND
TEND = tlim;

global UU
UU = U;




    %[FI] =computefluxint(ue,x,h,N,f,p, physics);
    FI = computeres(ue,x,h,N,f,p,physics,nSteps*k,gsp);




 global AD
 AD = computepseudo(N,x,h,r);


%global R
R = zeros(N+2,nSteps+1);
 tt=0;

for j = 1:nSteps+1
    

   R(:,j) =computeres(U(:,j),x,h,N,f,r,physics,tt,gsp);

tt = tt+k;
   

end




if(exist('goal','var') && strcmp(goal,'FI')==1)
    assert(strcmp(physics,'Poisson')==1)
for j = 1:nSteps+1
    R(:,j) = -FI;
end
end


Rm=max(abs(R(:,end)))



T=(0:1:nSteps)*k;
for j = 2:N+1
sp = spapi(6,T,R(j,:));
Rsp(j) = sp;
end


e = zeros(N+2,1);
ee =zeros(N+2,1);




global M
M = nSteps;


 AD = computepseudo(N,x,h,q);

 
 
T = 1;
s=1;

for j = 1:nSteps


E(:,j) = e;

    TT = k*(j-1);
    



% if( ((max(s)*k*inf<1e-15)||(TT>=tlim)) || (j >= nSteps))
if( ((TT>=tlim)) || (j >= nSteps))
    
[ee,s] = update('error',e,x,-R(:,j),k,h,N,q,tord,physics,TT,Rsp);



e = ee;
   max(s)
    TT
    T = (1:1:j)*k;
j
    break
end

s=0;

[ee,s] = update('error',e,x,-R(:,j),k,h,N,q,tord,physics,TT,Rsp);


e = ee;

T = (1:1:j)*k;


if(mod(j,100)==0)
max(s)
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



save('t','exacterr','ee','x')


end


clear global
end
