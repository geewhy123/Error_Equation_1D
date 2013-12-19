clear all
close all
N =80;
rng(1234);
h0 =1/N;
k = .0006*(40/N)^2;%0.00006;
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
    %f(i) = (1/h(i))*( log(1+(x(i)+h(i)/2))+((x(i)+h(i)/2)-1)/((x(i)+h(i)/2)+1) - (log(1+(x(i)-h(i)/2))+((x(i)-h(i)/2)-1)/((x(i)-h(i)/2)+1)) );%
    f(i)=(1/h(i))*(-2*(x(i)+h(i)/2)*exp(-(x(i)+h(i)/2)^2)+2*(x(i-1)+h(i-1)/2)*exp(-(x(i-1)+h(i-1)/2)^2));
    
    end
    %u(i) = (1/h(i))*((x(i)+h(i)/2)^4-((x(i)-h(i)/2)^4))/4;%exp(-(x(i)-0.5)^2);
    %u(i) = (1/h(i))*((x(i)+h(i)/2)^5-((x(i)-h(i)/2)^5))/5;%exp(-(x(i)-0.5)^2);
    %u(i) = (1/h(i))*((0.5*(1+(x(i)+h(i)/2))^2*log(1+(x(i)+h(i)/2))+7/4+1.5*(x(i)+h(i)/2)-0.25*(x(i)+h(i)/2)^2-2*(1+(x(i)+h(i)/2))*log(1+(x(i)+h(i)/2)))-((0.5*(1+(x(i)-h(i)/2))^2*log(1+(x(i)-h(i)/2))+7/4+1.5*(x(i)-h(i)/2)-0.25*(x(i)-h(i)/2)^2-2*(1+(x(i)-h(i)/2))*log(1+(x(i)-h(i)/2)))));
    u(i)=(1/h(i))*((pi^(1/2)*erf(x(i)+h(i)/2 ))/2-(pi^(1/2)*erf(x(i)-h(i)/2 ))/2);
    
end

f(1) = NaN;
f(N+2) = NaN;
 u(1) = NaN;
 u(N+2) = NaN;

 ue = u;
%u(:)=0;

%f(1)=2*1;
%f(N+2)=2*exp(-1);

%plot(x,u,'o')

%hold on

uu = zeros(N+2,1);
uu(1) = NaN;
uu(N+2) = NaN;
global AD
AD = computepseudo(N,x,h);
d =1;
figure
hold on
tic
olderr = 1;
%err= zeros(j,
T = 1;
for j = 1:50000
    
[err(j),Z]=unstructuredrecon1(u,x,h,N,1,exp(-1));
%if(j > 10 && ((abs(err(j)-olderr)/olderr < 1e-8) || abs(err(j)-olderr)<1e-15))
if (d*k < 1e-15)
    T = (1:1:j)*k;
    break
end

d=0;

olderr = err(j);
for i= 2:N+1
    y = Z(:,i);
    yr = Z(:,i+1);
    yl = Z(:,i-1);
    %need averaging flux

    
upr1 = y(2);
upr2 = yr(2);
ur = yr(1)+yr(2)*(-h(i)/2);
ul = y(1) + y(2)*(h(i)/2);
upr = (upr1+upr2)/2 + (0.2/h(i))*(ur(1)-ul(1));
upl1 = y(2);
upl2 = yl(2);

ur = y(1)+y(2)*(-h(i)/2);
ul = yl(1) + yl(2)*(h(i)/2);
upl = (upl1+upl2)/2 + (0.2/h(i))*(ur(1)-ul(1));


% upr1 = y(2);
% upr2 = yr(2);
% upr = (upr1+upr2)/2;% + 0.2*(yr(1)-y(1));
% upl1 = y(2);
% upl2 = yl(2);
% upl = (upl1+upl2)/2;% + 0.2*(y(1)-yl(1));

if i==N/2+1
   fhalf = upr; 
end


if i==2
    upl = upl1;
end
if i==N+1
   upr = upr1; 
end

uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
d = max(d,abs((upr-upl)/h(i)-f(i)));
end

u = uu;

T = (1:1:j)*k;
d
Z;
%plot(x,u,'*-')

end

shalf=exp(-0.5^2) - u(N/2+1)
fhalf=-exp(-0.5^2) - fhalf



cverr = sqrt(sum((ue(2:N+1)-u(2:N+1)).^2))/N

max(abs(ue-u))
%norm(ue(2:N+1)-u(2:N+1),2)
toc
lu20 = u;
lT20 = T;
lerr20 = err;
lx20 = x;
save('test','lu20','lT20','lerr20','lx20')
