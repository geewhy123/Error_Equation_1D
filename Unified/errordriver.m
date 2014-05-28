
function [errerr2,x,cverr2,exacterr,ee  ] = errordriver( N,p,q,r ,unif,BCLeft,valLeft,BCRight,valRight,tlim,tord,physics,goal)
%DRIVER Summary of this function goes here
%   Detailed explanation goes here

close all

if(p>0)

    
    
    rng(1234);

%  g = randi(1000000);
%  977219
% rng(g)
%   rng(972219);  
w = 0;
h0 = 1/N;

CFL = 0.2;
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
%    X(i) = X(i) + 0.1*randn*h0;
   X(i) = X(i) + unif*(-1+rand*(2))*h0/3;%0.001*sin(2*pi*X(i));%
   end
end

x = zeros(N+2,1);
for i = 2:N+1
    x(i) = (X(i-1)+X(i))/2;
end

x(1) = 0-(1-x(N+1));%-x(2);
x(N+2) = 1+x(2);%1+(1-x(N+1));
x

h = zeros(N+2,1);
for i = 2:N+1
   h(i) = X(i)-X(i-1); 
end

h(1) = h(N+1);
h(N+2) = h(2);


% global xx 
% xx = X;
global dir
dir = NaN*ones(N+2,1);


problem = pde(N,p,q,r,BCLeft,valLeft,BCRight,valRight,tlim,tord,physics,goal,x,h,k);
problem.initializeexact();

 
u0 = problem.initialSolution;
ue = problem.exactSolution;
% plot(x,u0,x,ue)
% assert(0)
problem.computemoments();

%u = ue;
u=u0;
%uu = zeros(N+2,1);

problem.exactSolution


%   if((strcmp(physics,'Poisson')==1 && strcmp(goal,'SS')==1 && problem.bcLeftType == 'D' && problem.bcRightType == 'D' )||(strcmp(physics,'Advection')==1 && strcmp(goal,'SS')==1))
 if(strcmp(goal,'SS')==1 && 0)
    fprintf('solving by Jacobian');
    [errerr2,x,cverr2,exacterr,ee  ]= problem.solvebyjacobian();
    
   
    return;
end



% global AD
% AD = computepseudo(N,x,h,p,BCRight,BCLeft);    
problem.computeprimalpseudo();

% problem.primalPI




 Z = problem.unstructuredrecon(ue,problem.pOrder,'solution');
 
er = problem.reconplot(Z)

f = problem.source;
% figure
% plot(x,ue)
% error('1')


% % % % [tau]=reconfluxsoln(Z,f,h,N,p,physics,tlim,problem)
[tau]=problem.computefluxintegral(Z,'solution')%reconfluxsoln(Z,f,h,N,p,physics,tlim,problem)


% tau

% error('1')
% max(abs(tau))
% % assert(0)
%  error('1')

%  computejacobiananalytic(p,h,N);
%   error('1')

% % % J = problem.computefluxjacobian(ue,'solution');%,x,h,N,p);
% % % Je = J
% % % % error('1')
% % % % % J
% % % % error('1')
% % % %   plot(x,J,x,f)
% % % 
% % % 
% % % % max(abs(J(2:N+1)-f(2:N+1))) 
% % % 
% % % % R = computeres(u,x,h,N,f,p,physics,tlim,NaN)
% % % [Z] = problem.unstructuredrecon(ue,p,'solution');%ue,x,h,N,NaN,NaN,p);
% % % 
% % % %  [er]=reconplot(x,h,N,p,Z);
% % % f = problem.source;
% % %  [tau]=reconfluxsoln(Z,f,h,N,p,physics,tlim,problem)
% % % 
% % %  max(abs(tau))
% % %   
% % %  del = ones(N,1);
% % %  t=0;
% % %  u0 = ue;
% % %  count = 0;
% % %  while(max(abs(del)) > 1e-10 && 0)
% % %      count = count +1;
% % %      dt = .01;
% % %      
% % % K = J(2:N+1,2:N+1)+eye(N)/dt;
% % % % K
% % % % error('1')
% % %  K = (K+K')/2;
% % % 
% % % [Z] = problem.unstructuredrecon(u0,p,'solution');%u0,x,h,N,NaN,NaN,p);
% % % 
% % % %  [er]=reconplot(x,h,N,p,Z);
% % %  [R]=reconfluxsoln(Z,f,h,N,p,physics,t,problem)
% % %     del = K\R(2:N+1);
% % %      uu = u0(2:N+1) + del*dt;
% % %      u0 = NaN*ones(N+2,1);
% % %      u0(2:N+1) = uu;
% % %      t = t+dt;
% % %  end
% % %  u0
% % %  max(abs(u0-ue))
% % %  u0-ue
% % % 
% % %  v=Je(2:N+1,2:N+1)\tau(2:N+1)
% % %  figure
% % %  plot(x,u0-ue,x(2:N+1),v)
% % %  figure
%   error('1')

 

d=1;

T = 1;
for j = 1:100000
        U(:,j) = u;
tt = k*(j-1);
    


if((max(d)*k<1e-15)||(tt>=tlim))
     
%  [uu,d] = update('solution',u,x,f,k,h,N,p,tord,physics,NaN,NaN);
u = uu;
    max(d)
    tt
    T = (1:1:j-1)*k;
%  U(:,j+1) = u;
 nSteps = j-1;
 
    break
end

d=0;


[uu,d] = problem.updatesolution(u);

% [uu,d] = update('solution',u,x,problem.source,k,h,N,p,tord,physics,NaN,NaN,problem);

% if(j==20)
%     uu
% error('1')
% end
% [uu,d] = update('solution',u,x,f,k,h,N,p,tord,physics,NaN,NaN);


u = uu;

% T = (1:1:j)*k;

if(mod(j,100)==0)
    max(d)
end

end

u

cverr1 = sum(abs(ue(2:N+1)-u(2:N+1)))/N
cverr2 = sqrt(sum((ue(2:N+1)-u(2:N+1)).^2)/N)
cverrinf=max(abs(ue-u))
size(U)


u(1) = NaN;
u(N+2) = NaN;
plot(x,u,'*',x,ue,'o')
ylabel('u')
figure
plot(x,ue-u,'x')
ylabel('ue-u')


end

T(end)
% 
% if( T(end) > 0.997)%0.9166)
% g
% error('1')
% else
% errerr2=NaN;
% ee = NaN;
% exacterr = NaN;
%     return
% end

assert((abs(T(end)-tlim)/tlim < 1e-4) || (strcmp(physics,'Poisson')==1 && tlim/T(end) > 2 ) ) 

tlim = T(end);




global dUdt
dUdt = diffU(U,k);

%  U
%  error('1')
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




% 
% t = (0:1:nSteps)*k;
% Utexact = zeros(size(dUdt));
% Utexact(1,:) = NaN;
% Utexact(N+2,:) = NaN;
% Uexact = Utexact;
% for i = 2:N+1
%     for j = 1:nSteps+1
%           xl = x(i)-h(i)/2;
%     xr = x(i)+h(i)/2;
%         Uexact(i,j) = -(1/(2*pi))*(1/h(i))*(cos(2*pi*(xr+t(j)))-cos(2*pi*(xl+t(j))));
%         Utexact(i,j) =           (1/h(i))*(sin(2*pi*(xr+t(j)))-sin(2*pi*(xl+t(j))));
%     end
% end
% 
% dUdt = diffU(Uexact,k);
% 
% max(max(abs(dUdt-Utexact)))
% error('1')


nSteps
% 
% if(nSteps < 20)
%    error('1') 
% end

dir

if(q>0 && r > 0)
    
    clearvars -except u N p q r unif FI bta f cverr2 v k ue u0 tlim tord uo physics uder nSteps gsp U h x goal dUdt X problem Je tau
    
    
    
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

% %     FI = computeres(ue,x,h,N,f,p,physics,nSteps*k,gsp);
% % %  FI = problem.computeres(ue,nSteps*k,p);
 
 
%     FI
% error('1')


%  global AD
  problem.computerespseudo();%N,x,h,r);
% problem.resPI
% error('1')

%global R
R = zeros(N+2,nSteps+1);
 tt=0;

for j = 1:nSteps+1
    

%    R(:,j) =computeres(U(:,j),x,h,N,f,r,physics,tt,gsp);
R(:,j) = problem.computeres(U(:,j),tt,r);


tt = tt+k;
   

end


if(exist('goal','var') && strcmp(goal,'FI')==1)
    assert(strcmp(physics,'Poisson')==1)
for j = 1:nSteps+1
    R(:,j) = -FI;
end
end



problem.residual = R;

[-R(:,end) tau]
max(abs(-R(:,end)-tau))
plot(x,-R(:,end),x,tau)
R
% error('1')






Rm=max(abs(R(:,end)))



 sqrt(sum((R(2:N+1,end)).^2)/N)
%  error('1');
% 
%     cverr2 = Rm;
%     errerr2 = Rm;
%     ee = NaN;
%     exacterr = NaN;
%     return


plot(x,R(:,end))
max(abs(R(:,end)))
%  error('2')

% clear R
%  Ro = R;
% R
% error('1')
%     load('RA.mat')
%  Ro(:,1) = R(:,1);
%  R = Ro;


% dUdt 
% 
% t = (0:1:nSteps)*k;
% Utexact = zeros(size(dUdt));
% Utexact(1,:) = NaN;
% Utexact(N+2,:) = NaN;
% Uexact = Utexact;
% for i = 2:N+1
%     for j = 1:nSteps+1
%           xl = x(i)-h(i)/2;
%     xr = x(i)+h(i)/2;
%         Uexact(i,j) = -(1/(2*pi))*(1/h(i))*(cos(2*pi*(xr+t(j)))-cos(2*pi*(xl+t(j))));
%         Utexact(i,j) =           (1/h(i))*(sin(2*pi*(xr+t(j)))-sin(2*pi*(xl+t(j))));
%     end
% end
% 
% Utexact
% plot(x,dUdt(:,end)-Utexact(:,end))
% max(abs(dUdt(:,end)-Utexact(:,end)))
% % plot(x,Uexact(:,end))
% 
%  error('1')



T=(0:1:nSteps)*k;

for j = 2:N+1
sp = spapi(6,T,R(j,:));
Rsp(j) = sp;
end


e = zeros(N+2,1);
ee =zeros(N+2,1);


problem.Rsp =Rsp;

global M
M = nSteps;


% % %  AD = computepseudo(N,x,h,q);
problem.computeerrorpseudo();



%new
% % % Je = problem.computefluxjacobian(ue,'error');
% % % % problem.errorRM
% % % % error('1')
% % % 
% % % [Z] = problem.unstructuredrecon(ue-u,q,'error');%ue,x,h,N,NaN,NaN,p);
% % % f = -R(:,end);
% % % 
% % % % error('1')
% % %  [tauE]=reconfluxsoln(Z,f,h,N,q,physics,tlim,problem)
% % % %  Je
% % % %  error('1')
% % % w = Je(2:N+1,2:N+1)\tauE(2:N+1);
% % % max(abs(w))
% % % % error('1')
% % % figure
% % % plot(x(2:N+1),w)
% % % figure
% % % % error('2')

%new






% size(U)
% size(R)

%  
T = 1;
s=1;

%  nSteps = nSteps-1;
% tlim = tlim-k;


for j = 1:nSteps+1


E(:,j) = e;

    TT = k*(j-1);
    



% if( ((max(s)*k*inf<1e-15)||(TT>=tlim)) || (j >= nSteps))
if( ((TT>=tlim)) || (j >= nSteps+1))
%     TT
%     tlim
 
%  [ee,s] = update('error',e,x,-R(:,j),k,h,N,q,tord,physics,TT,Rsp);



e = ee;
   max(s)
    TT
    T = (1:1:j-1)*k;
j
nSteps
% R(:,end-2:end)
    break
end

s=0;

% % % [ee,s] = update('error',e,x,-R(:,j),k,h,N,q,tord,physics,TT,Rsp);

[ee,s] = problem.updateerror(e,TT,j);


e = ee;

% T = (1:1:j-1)*k;


if(mod(j,100)==0)
max(s)
end
end


exacterr = ue-u;

% norm(u(2:N+1))

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
ylabel('\epsilon - \epsilon_h')


save('t','exacterr','ee','x')


else 
    errerr2 = NaN;
    ee = NaN;
    exacterr = NaN;
end


exacterr-ee

clear global
end
