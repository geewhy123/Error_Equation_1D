function [ output_args ] = unstructuredrecon1D( N,p )
%UNSTRUCTUREDRECON1 Summary of this function goes here
%   Detailed explanation goes here

close all
%N = 20;
%p=4;
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
   ;
    else
      xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;
    %f(i)= (1/h(i))*(1/(12*pi))*(3*cos(2*pi*xr)-cos(6*pi*xr)-3*cos(2*pi*xl)+cos(6*pi*xl));
    %f(i) = (1/h(i))*(-2*pi);%*(exp(cos(2*pi*xr))*   (2*(cos(pi*xr)^2-    4*(cos(pi*xr))^4+1)));% - (exp(cos(2*pi*xl))*(2*(cos(pi*xl))^2-4*(cos(pi*xl))^4+1));
    %f(i)=(1/h(i))*((-1/pi)*sin(2*pi*xr)+(1/(pi^3))*(sin(4*pi*xr))+(1/pi)*(sin(2*pi*xl))-(1/(pi^3))*(sin(4*pi*xl)));
  %f(i) = (1/h(i))*(-4*pi^2)*( (exp(1)^3*sin(2*pi*xr)+1)/(sin(2*pi*xr)+exp(1)^3)^2 - (exp(1)^3*sin(2*pi*xl)+1)/(sin(2*pi*xl)+exp(1)^3)^2);
   f(i) = (1/h(i))*(-2*pi)*(sin(2*pi*xr)-sin(2*pi*xl));
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
  %u(i) = (1/h(i))*(log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl)));
    u(i) = (1/h(i))*(1/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl)); 
    
    up(i) = (1/h(i))*(cos(2*pi*xr)-cos(2*pi*xl));
   
end

 u(1) = NaN;
 u(N+2) = NaN;

u
%f(1)=2*1;
%f(N+2)=2*exp(-1);



%plot(x,u,'o')

hold on 

global AD
AD = computepseudo(N,x,h,p);
%uu=zeros(N+2,1);
%for t = 1:10
%k = 0.1;
%u=uu;

err = 0;
M = 10000;
for i = 2:N+1

switch p
    case 2
        switch i
            case 2
                y= recon1(x(i),h(i),u(i),x(N+1),h(N+1),u(N+1),x(i+1),h(i+1),u(i+1),x(N),h(N),u(N),x(i+2),h(i+2),u(i+2),i)
            case 3
                y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(N+1),h(N+1),u(N+1),x(i+2),h(i+2),u(i+2),i)
            case N
                y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(2),h(2),u(2),i);
            case N+1
                y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(2),h(2),u(2),x(i-2),h(i-2),u(i-2),x(3),h(3),u(3),i);   
            otherwise
                y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),i);
        end
   
        xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
   ue = cos(2*pi*xx);%2*pi*cos(2*pi*xx)./(sin(2*pi*xx)+exp(1)^3);%exp(-(xx-0.5).^2);%xx.^3;%;exp(-(xx-0.5).^2);
   %ue = xx.^3;
   %ue = xx.^4;
   yy = y(1)+y(2)*(xx-x(i));%+y(3)*(xx-x(i)).^2;%+y(4)*(xx-x(i)).^3;%+y(5)*(xx-x(i)).^4;
   
    case 3
        switch i
            case 2
                y= recon2(x(i),h(i),u(i),x(N+1),h(N+1),u(N+1),x(i+1),h(i+1),u(i+1),x(N),h(N),u(N),x(i+2),h(i+2),u(i+2),i);
            case 3
                y= recon2(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(N+1),h(N+1),u(N+1),x(i+2),h(i+2),u(i+2),i);
            case N
                y= recon2(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(2),h(2),u(2),i);
            case N+1
                y= recon2(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(2),h(2),u(2),x(i-2),h(i-2),u(i-2),x(3),h(3),u(3),i);   
            otherwise
                y= recon2(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),i);
        end
          Z(:,i) = y;
   
        xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
   ue = cos(2*pi*xx);%2*pi*cos(2*pi*xx)./(sin(2*pi*xx)+exp(1)^3);%exp(-(xx-0.5).^2);%xx.^3;%;exp(-(xx-0.5).^2);
   %ue = xx.^3;
   %ue = xx.^4;
   yy = y(1)+y(2)*(xx-x(i))+y(3)*(xx-x(i)).^2;%+y(4)*(xx-x(i)).^3;%+y(5)*(xx-x(i)).^4;

%   if(i==2)
 %      h(i)^(2)/(3*2^2)
   ub(i) = y(1)+y(3)*h(i)^(2)/(3*2^2);
   upb(i) = y(2);
   uppb(i) = y(3);
   %end
%    if(i==2)
%        h(i)^(3)/(3*2^2)
%       y
%       ub(i) 
%    end
   
    case 4
        switch i
            case 2
                y= recon3(x(i),h(i),u(i),x(N+1),h(N+1),u(N+1),x(i+1),h(i+1),u(i+1),x(N),h(N),u(N),x(i+2),h(i+2),u(i+2),i)
            case 3
                y= recon3(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(N+1),h(N+1),u(N+1),x(i+2),h(i+2),u(i+2),i)
            case N
                y= recon3(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(2),h(2),u(2),i);
            case N+1
                y= recon3(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(2),h(2),u(2),x(i-2),h(i-2),u(i-2),x(3),h(3),u(3),i);   
            otherwise
                y= recon3(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),i);
        end
      Z(:,i) = y;
        xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
   ue = cos(2*pi*xx);%2*pi*cos(2*pi*xx)./(sin(2*pi*xx)+exp(1)^3);%exp(-(xx-0.5).^2);%xx.^3;%;exp(-(xx-0.5).^2);
   %ue = xx.^3;
   %ue = xx.^4;
   yy = y(1)+y(2)*(xx-x(i))+y(3)*(xx-x(i)).^2+y(4)*(xx-x(i)).^3;%+y(5)*(xx-x(i)).^4;

   
   %ub(i) = y(1)+y(3)*h(i)^(4)/(5*2^4);

  case 5
        switch i
            case 2
                y= recon4(x(i),h(i),u(i),x(N+1),h(N+1),u(N+1),x(i+1),h(i+1),u(i+1),x(N),h(N),u(N),x(i+2),h(i+2),u(i+2),i);
            case 3
                y= recon4(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(N+1),h(N+1),u(N+1),x(i+2),h(i+2),u(i+2),i);
            case N
                y= recon4(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(2),h(2),u(2),i);
            case N+1
                y= recon4(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(2),h(2),u(2),x(i-2),h(i-2),u(i-2),x(3),h(3),u(3),i);   
            otherwise
                y= recon4(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),i);
        end
        
      Z(:,i)=y;
   
        xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
   ue = cos(2*pi*xx);%2*pi*cos(2*pi*xx)./(sin(2*pi*xx)+exp(1)^3);%exp(-(xx-0.5).^2);%xx.^3;%;exp(-(xx-0.5).^2);
   %ue = xx.^3;
   %ue = xx.^4;
   yy = y(1)+y(2)*(xx-x(i))+y(3)*(xx-x(i)).^2+y(4)*(xx-x(i)).^3+y(5)*(xx-x(i)).^4;
  
   case 6
        switch i
        
    case 2
 y= recon5(x(i),h(i),u(i),x(N+1),h(N+1),u(N+1),x(i+1),h(i+1),u(i+1),x(N),h(N),u(N),x(i+2),h(i+2),u(i+2),x(N-1),h(N-1),u(N-1),x(i+3),h(i+3),u(i+3),i);
    case 3
  y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(N+1),h(N+1),u(N+1),x(i+2),h(i+2),u(i+2),x(N),h(N),u(N),x(i+3),h(i+3),u(i+3),i);
    case 4
         y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),x(N+1),h(N+1),u(N+1),x(i+3),h(i+3),u(i+3),i);
    case N-1
         y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),x(i-3),h(i-3),u(i-3),x(2),h(2),u(2),i);
    case N
 y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(2),h(2),u(2),x(i-3),h(i-3),u(i-3),x(3),h(3),u(3),i);
    case N+1
   y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(2),h(2),u(2),x(i-2),h(i-2),u(i-2),x(3),h(3),u(3),x(i-3),h(i-3),u(i-3),x(4),h(4),u(4),i);
    otherwise
  y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),x(i-3),h(i-3),u(i-3),x(i+3),h(i+3),u(i+3),i);
        end
      
   
        xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
   ue = cos(2*pi*xx);%2*pi*cos(2*pi*xx)./(sin(2*pi*xx)+exp(1)^3);%exp(-(xx-0.5).^2);%xx.^3;%;exp(-(xx-0.5).^2);
   %ue = xx.^3;
   %ue = xx.^4;
   yy = y(1)+y(2)*(xx-x(i))+y(3)*(xx-x(i)).^2+y(4)*(xx-x(i)).^3+y(5)*(xx-x(i)).^4;
   
   
end
   plot(xx,yy)
err = max(err,max(abs(ue-yy)));
end


%Z
%Z
u

ub(N+2) = NaN;
ub-u
upb(N+2) = NaN;
upb(1) = NaN;
d1 = max(abs(upb-up))


uppb(1) = NaN;
uppb(N+2) = NaN;
d2 = max(abs(uppb'-f))

 % uberr=max(ub-u) 
 
% uperr = max(
   err
%end
  

end

