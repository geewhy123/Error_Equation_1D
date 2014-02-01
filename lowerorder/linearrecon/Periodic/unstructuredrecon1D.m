clear all
close all
N = 1280;
p=2;
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
  f(i) = (1/h(i))*(-4*pi^2)*( (exp(1)^3*sin(2*pi*xr)+1)/(sin(2*pi*xr)+exp(1)^3)^2 - (exp(1)^3*sin(2*pi*xl)+1)/(sin(2*pi*xl)+exp(1)^3)^2);
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
  u(i) = (1/h(i))*(log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl)));
    
end

 u(1) = NaN;
 u(N+2) = NaN;


f(1)=2*1;
f(N+2)=2*exp(-1);

%plot(x,u,'o')

hold on 

global AD
AD = computepseudo(N,x,h,p);
%uu=zeros(N+2,1);
%for t = 1:10
%k = 0.1;
%u=uu;

% %second order recon
% for i = 2:N+1
% %i = 5;
% wi1 = 1/abs(x(i-1)-x(i));
% wi2 = 1/abs(x(i+1)-x(i));
% xi = x(i);
% x1 = x(i-1);
% x2 = x(i+1);
% syms z
% xb1 = (1/h(i-1))*int(z-x(i-1),z,x(i-1)-h(i-1)/2,x(i-1)+h(i-1)/2);
% xb2 = (1/h(i+1))*int(z-x(i+1),z,x(i+1)-h(i+1)/2,x(i+1)+h(i+1)/2);
% xbi = (1/h(i))*int(z-x(i),z,x(i)-h(i)/2,x(i)+h(i)/2);
% ub1 = u(i-1);
% ub2 = u(i+1);
% ubi = u(i);
% A = double([wi1*(xb1-xbi+x1-xi); wi2*(xb2-xbi+x2-xi)]);
% b = [wi1*(ub1-ubi); wi2*(ub2-ubi)];
% 
% y(2) = A\b%(A'*A)\(A'*b)
% y(1) = ubi%ubi-xbi*y(2)
% 
% 
%    xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,100);
%    yy = y(1)+y(2)*(xx-x(i));
%    plot(xx,yy)
% end
error = 0;
M = 10000;
%4th order recon
for i = 2:N+1
%i = 5;
% wi1 = 1/abs(x(i-1)-x(i));
% wi2 = 1/abs(x(i+1)-x(i));
% wi3 = 1/abs(x(i-2)-x(i));
% wi4 = 1/abs(x(i+2)-x(i));
% xi = x(i);
% x1 = x(i-1);
% x2 = x(i+1);
% x3 = x(i-2);
% x4 = x(i+2);
% syms z
% xb1 = (1/h(i-1))*int(z-x(i-1),z,x(i-1)-h(i-1)/2,x(i-1)+h(i-1)/2);
% xb2 = (1/h(i+1))*int(z-x(i+1),z,x(i+1)-h(i+1)/2,x(i+1)+h(i+1)/2);
% xb3 = (1/h(i-2))*int(z-x(i-2),z,x(i-2)-h(i-2)/2,x(i-2)+h(i-2)/2);
% xb4 = (1/h(i+2))*int(z-x(i+2),z,x(i+2)-h(i+2)/2,x(i+2)+h(i+2)/2);
% xbi = (1/h(i))*int(z-x(i),z,x(i)-h(i)/2,x(i)+h(i)/2);
% 
% x2b1 = (1/h(i-1))*int((z-x(i-1))^2,z,x(i-1)-h(i-1)/2,x(i-1)+h(i-1)/2);
% x2b2 = (1/h(i+1))*int((z-x(i+1))^2,z,x(i+1)-h(i+1)/2,x(i+1)+h(i+1)/2);
% x2b3 = (1/h(i-2))*int((z-x(i-2))^2,z,x(i-2)-h(i-2)/2,x(i-2)+h(i-2)/2);
% x2b4 = (1/h(i+2))*int((z-x(i+2))^2,z,x(i+2)-h(i+2)/2,x(i+2)+h(i+2)/2);
% x2bi = (1/h(i))*int((z-x(i))^2,z,x(i)-h(i)/2,x(i)+h(i)/2);
% 
% x3b1 = (1/h(i-1))*int((z-x(i-1))^3,z,x(i-1)-h(i-1)/2,x(i-1)+h(i-1)/2);
% x3b2 = (1/h(i+1))*int((z-x(i+1))^3,z,x(i+1)-h(i+1)/2,x(i+1)+h(i+1)/2);
% x3b3 = (1/h(i-2))*int((z-x(i-2))^3,z,x(i-2)-h(i-2)/2,x(i-2)+h(i-2)/2);
% x3b4 = (1/h(i+2))*int((z-x(i+2))^3,z,x(i+2)-h(i+2)/2,x(i+2)+h(i+2)/2);
% x3bi = (1/h(i))*int((z-x(i))^3,z,x(i)-h(i)/2,x(i)+h(i)/2);
% ub1 = u(i-1);
% ub2 = u(i+1);
% ub3 = u(i-2);
% ub4 = u(i+2);
% ubi = u(i);
% A = double([wi1*(xb1-xbi+x1-xi) wi1*(x2b1+2*(x1-xi)*xb1+(x1-xi)^2-xbi) wi1*(x3b1+3*(x1-xi)*x2b1+3*(x1-xi)^2*xb1+(x1-xi)^3-xbi); 
%             wi2*(xb2-xbi+x2-xi) wi2*(x2b2+2*(x2-xi)*xb2+(x2-xi)^2-xbi) wi2*(x3b2+3*(x2-xi)*x2b2+3*(x2-xi)^2*xb2+(x2-xi)^3-xbi); 
%             wi3*(xb3-xbi+x3-xi) wi3*(x2b3+2*(x3-xi)*xb3+(x3-xi)^2-xbi) wi3*(x3b3+3*(x3-xi)*x2b3+3*(x3-xi)^2*xb3+(x3-xi)^3-xbi);
%             wi4*(xb4-xbi+x4-xi) wi4*(x2b4+2*(x4-xi)*xb4+(x4-xi)^2-xbi) wi4*(x3b4+3*(x4-xi)*x2b4+3*(x4-xi)^2*xb4+(x4-xi)^3-xbi)]);
% b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ];
% 
% y(2:4) = A\b;%(A'*A)\(A'*b)
% y(1) = ubi-xbi*y(2)-x2bi*y(3)-x3bi*y(4);%ubi-xbi*y(2)
% %q = y(1)-ubi


switch i
        case 2
 y= recon1(x(i),h(i),u(i),x(N+1),h(N+1),u(N+1),x(i+1),h(i+1),u(i+1),x(N),h(N),u(N),x(i+2),h(i+2),u(i+2),i);
    case 3
  y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(N+1),h(N+1),u(N+1),x(i+2),h(i+2),u(i+2),i);
    case N
 y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(2),h(2),u(2),i);
    case N+1
   y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(2),h(2),u(2),x(i-2),h(i-2),u(i-2),x(3),h(3),u(3),i);
    
    otherwise
 y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),i);
end
   
   
   xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
   ue = 2*pi*cos(2*pi*xx)./(sin(2*pi*xx)+exp(1)^3);%exp(-(xx-0.5).^2);%xx.^3;%;exp(-(xx-0.5).^2);
   %ue = xx.^3;
   %ue = xx.^4;
   yy = y(1)+y(2)*(xx-x(i));%+y(3)*(xx-x(i)).^2+y(4)*(xx-x(i)).^3;
   plot(xx,yy)
   %ylim([-1e0 1e0])
  
   
   error = max(error,max(abs(ue-yy)));
   %upr = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
   %upl = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
   %uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
end


% % i = 2
% i=2;
% 
% y= recon3(x(i),h(i),u(i),x(i+1),h(i+1),u(i+1),x(i+2),h(i+2),u(i+2),x(i+3),h(i+3),u(i+3),x(i+4),h(i+4),u(i+4),i);
% 
%    xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
%    yy = y(1)+y(2)*(xx-x(i))+y(3)*(xx-x(i)).^2+y(4)*(xx-x(i)).^3;
% 
%    plot(xx,yy)
%    %ue = xx.^3;
%    ue = exp(-(xx-0.5).^2);%xx.^3;
%    error = max(error,max(abs(ue-yy)));
%    %upr = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
%    %upl = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
%    %uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
%   
% 
%    
% 
%    
%    i = N+1;
%    y= recon3(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i-2),h(i-2),u(i-2),x(i-3),h(i-3),u(i-3),x(i-4),h(i-4),u(i-4),i);
% 
%    xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
%    yy = y(1)+y(2)*(xx-x(i))+y(3)*(xx-x(i)).^2+y(4)*(xx-x(i)).^3;
%    plot(xx,yy)
% %   ue = xx.^3;
%    ue = exp(-(xx-0.5).^2);%xx.^3;
%    error = max(error,max(abs(ue-yy)));
%    %upr = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
%    %upl = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
%    %uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
%    
%    
%    i=3;
% 
% y= recon3(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i+2),h(i+2),u(i+2),x(i+3),h(i+3),u(i+3),i);
% 
%    xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
%    yy = y(1)+y(2)*(xx-x(i))+y(3)*(xx-x(i)).^2+y(4)*(xx-x(i)).^3;
%    plot(xx,yy)
% %ue = xx.^3;
% ue = exp(-(xx-0.5).^2);%xx.^3;
%    error = max(error,max(abs(ue-yy)));
%    %upr = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
%    %upl = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
%    %uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
% 
%    
%    i = N;
%    y= recon3(x(i),h(i),u(i),x(i+1),h(i+1),u(i+1),x(i-1),h(i-1),u(i-1),x(i-2),h(i-2),u(i-2),x(i-3),h(i-3),u(i-3),i);
% 
%    xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
%    yy = y(1)+y(2)*(xx-x(i))+y(3)*(xx-x(i)).^2+y(4)*(xx-x(i)).^3;
%    plot(xx,yy)
% 
% %ue = xx.^3;
% ue = exp(-(xx-0.5).^2);%xx.^3;
%    error = max(error,max(abs(ue-yy)))
%   %  upr = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
%   % upl = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
%   % uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
   
   error
%end
  