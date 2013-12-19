clear all
close all
N = 10;
rng(1234);
h0 = 1/N;
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
for i = 1:N+2
    if i==1
       f(i) = 0; 
    else
    f(i) = (1/h(i))*(-2*(x(i)+h(i)/2)*exp(-(x(i)+h(i)/2)^2)+2*(x(i-1)+h(i-1)/2)*exp(-(x(i-1)+h(i-1)/2)^2));
    
    end
  %  u(i) = (1/h(i))*((x(i)+h(i)/2)^4-((x(i)-h(i)/2)^4))/4;%exp(-(x(i)-0.5)^2);
   %u(i) = (1/h(i))*((x(i)+h(i)/2)^5-((x(i)-h(i)/2)^5))/5;%exp(-(x(i)-0.5)^2);
%   u(i) = (1/h(i))*(sin(pi*(x(i)+h(i)/2))/pi - sin(pi*(x(i)-h(i)/2))/pi);
  % u(i) = (1/h(i))*((pi^(1/2)*erf(x(i)+h(i)/2 - 1/2))/2-(pi^(1/2)*erf(x(i)-h(i)/2 - 1/2))/2);
     u(i)=(1/h(i))*((pi^(1/2)*erf(x(i)+h(i)/2 ))/2-(pi^(1/2)*erf(x(i)-h(i)/2 ))/2);
end

 u(1) = NaN;
 u(N+2) = NaN;


f(1)=2*1;
f(N+2)=2*exp(-1);

%plot(x,u,'o')

hold on 

global AD
AD = computepseudo(N,x,h);

error = 0;
M = 10000;
%4th order recon
for i = 4:N-1


 y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),i);
   
   
   
   xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
   ue = exp(-xx.^2);%cos(pi*xx);%exp(-(xx-0.5).^2);%xx.^3;%;exp(-(xx-0.5).^2);
   %ue = xx.^3;
   %ue = xx.^4;
   yy = y(1)+y(2)*(xx-x(i));
   plot(xx,yy)
   %ylim([-1e0 1e0])
  
   
   error = max(error,max(abs(ue-yy)));
   %upr = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
   %upl = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
   %uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
end


% i = 2
i=2;

%%%y= recon1(x(i),h(i),u(i),x(i+1),h(i+1),u(i+1),x(i+2),h(i+2),u(i+2),x(i+3),h(i+3),u(i+3),x(i+4),h(i+4),u(i+4),i);
y= reconboundary1(x(i),h(i),u(i),x(i+1),h(i+1),u(i+1),x(i+2),h(i+2),u(i+2),x(i+3),h(i+3),u(i+3),x(i+4),h(i+4),u(i+4),1,'left');

xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
   yy = y(1)+y(2)*(xx-x(i));

   plot(xx,yy)
  % ue = xx.^4;
   ue = exp(-xx.^2);%cos(pi*xx);%exp(-(xx-0.5).^2);%xx.^3;
   error = max(error,max(abs(ue-yy)));
   %upr = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
   %upl = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
   %uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
  

   

   
   i = N+1;
   %%%y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i-2),h(i-2),u(i-2),x(i-3),h(i-3),u(i-3),x(i-4),h(i-4),u(i-4),i);
   y= reconboundary1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i-2),h(i-2),u(i-2),x(i-3),h(i-3),u(i-3),x(i-4),h(i-4),u(i-4),exp(-1),'right');

   xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
   yy = y(1)+y(2)*(xx-x(i));
   plot(xx,yy)
 % ue = xx.^4;
   ue = exp(-xx.^2);%cos(pi*xx);%exp(-(xx-0.5).^2);%xx.^3;
   error = max(error,max(abs(ue-yy)));
   %upr = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
   %upl = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
   %uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
   
   
   i=3;

y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i+2),h(i+2),u(i+2),x(i+3),h(i+3),u(i+3),i);

   xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
   yy = y(1)+y(2)*(xx-x(i));
   plot(xx,yy)
%ue = xx.^4;
ue = exp(-xx.^2);%cos(pi*xx);%exp(-(xx-0.5).^2);%xx.^3;
   error = max(error,max(abs(ue-yy)));
   %upr = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
   %upl = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
   %uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 

   
   i = N;
   y= recon1(x(i),h(i),u(i),x(i+1),h(i+1),u(i+1),x(i-1),h(i-1),u(i-1),x(i-2),h(i-2),u(i-2),x(i-3),h(i-3),u(i-3),i);

   xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
   yy = y(1)+y(2)*(xx-x(i));
   plot(xx,yy)

%ue = xx.^4;
ue =exp(-xx.^2);% cos(pi*xx);%exp(-(xx-0.5).^2);%xx.^3;
   error = max(error,max(abs(ue-yy)))
  %  upr = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
  % upl = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
  % uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
   
   
%end
  