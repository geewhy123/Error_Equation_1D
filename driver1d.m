clear all
close all
N = 20;
rng(1234);
h0 = 1/N;
k = .0008%0.00006;
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

u = zeros(N+2,1);

f = zeros(N+2,1);
for i = 1:N+2
    if i==1
       f(i) = 0; 
    else
    f(i) = (1/h(i))*(-2*(x(i)+h(i)/2)*exp(-(x(i)+h(i)/2)^2)+2*(x(i-1)+h(i-1)/2)*exp(-(x(i-1)+h(i-1)/2)^2));
    
    end
    %u(i) = (1/h(i))*((x(i)+h(i)/2)^4-((x(i)-h(i)/2)^4))/4;%exp(-(x(i)-0.5)^2);
    %u(i) = (1/h(i))*((x(i)+h(i)/2)^5-((x(i)-h(i)/2)^5))/5;%exp(-(x(i)-0.5)^2);
    u(i) = (1/h(i))*((pi^(1/2)*erf(x(i)+h(i)/2 ))/2-(pi^(1/2)*erf(x(i)-h(i)/2 ))/2);
    
end

 u(1) = NaN;
 u(N+2) = NaN;

 
%u(:)=0;

f(1)=2*1;
f(N+2)=2*exp(-1);

%plot(x,u,'o')

%hold on

uu = zeros(N+2,1);

global AD
AD = computepseudo(N,x,h);


tic
olderr = 1;
%err= zeros(j,
T = 1;
for j = 1:100000
    
[err(j),Z]=unstructuredrecon3(u,x,h,N,1,exp(-1));
if(j > 10 && ((abs(err(j)-olderr)/olderr < 1e-8) || abs(err(j)-olderr)<1e-15))
    T = (1:1:j)*k;
    break
end
olderr = err(j);
for i= 2:N+1
    y = Z(:,i);
    yr = Z(:,i+1);
    yl = Z(:,i-1);
    %need averaging flux
upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2;
upr = (upr1+upr2)/2;
upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2;
upl = (upl1+upl2)/2;

if i==2
    upl = upl1;
end
if i==N+1
   upr = upr1; 
end

uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
end

u = uu;

T = (1:1:j)*k;
end
toc
nu10 = u;
nT10 = T;
nerr10 = err;
nx10 = x;
save('test10','nu10','nT10','nerr10','nx10')
