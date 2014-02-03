function [error, Z] = unstructuredrecon4( u,x,h,N,u0,u1 )
%UNSTRUCTUREDRECON3 Summary of this function goes here
%   Detailed explanation goes here
Z = zeros(5,N+2);
error = 0;
M = 1000;
%4th order recon
for i = 2:N+1%4:N-1

switch i
    case 2
 Y= recon4(x(i),h(i),u(i),x(N+1),h(N+1),u(N+1),x(i+1),h(i+1),u(i+1),x(N),h(N),u(N),x(i+2),h(i+2),u(i+2),i);
    case 3
  Y= recon4(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(N+1),h(N+1),u(N+1),x(i+2),h(i+2),u(i+2),i);
    case N
 Y= recon4(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(2),h(2),u(2),i);
    case N+1
   Y= recon4(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(2),h(2),u(2),x(i-2),h(i-2),u(i-2),x(3),h(3),u(3),i);
    otherwise
  Y= recon4(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),i);
 
   
   
   
%   xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
%   ue = exp(-(xx).^2);%xx.^3;%;exp(-(xx-0.5).^2);
   %ue = xx.^3;
   %ue = xx.^4;
%   yy = y(1)+y(2)*(xx-x(i))+y(3)*(xx-x(i)).^2+y(4)*(xx-x(i)).^3+y(5)*(xx-x(i)).^4;
  % plot(xx,yy)
   %ylim([-1e0 1e0])
  
   
 %  error = max(error,max(abs(ue-yy)));
   %upr = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
   %upl = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
   %uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
  
end

 y=Y;
   Z(:,i) = y;

% % i = 2
% i=2;
% 
% y= reconboundary4(x(i),h(i),u(i),x(i+1),h(i+1),u(i+1),x(i+2),h(i+2),u(i+2),x(i+3),h(i+3),u(i+3),x(i+4),h(i+4),u(i+4),u0,'left');
% 
%   % xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
%   % yy = y(1)+y(2)*(xx-x(i))+y(3)*(xx-x(i)).^2+y(4)*(xx-x(i)).^3+y(5)*(xx-x(i)).^4;
% %ue = exp(-(xx).^2);%xx.^3;
%  %  plot(xx,yy)
%    %ue = xx.^3;
%    
%  %  error = max(error,max(abs(ue-yy)));
%    %upr = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
%    %upl = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
%    %uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
%   Z(:,i)= y;
% 
%    
% 
%    
%    i = N+1;
%    y= reconboundary4(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i-2),h(i-2),u(i-2),x(i-3),h(i-3),u(i-3),x(i-4),h(i-4),u(i-4),u1,'right');
% 
%   % xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
%   % yy = y(1)+y(2)*(xx-x(i))+y(3)*(xx-x(i)).^2+y(4)*(xx-x(i)).^3+y(5)*(xx-x(i)).^4;
%   %  ue = exp(-(xx).^2);%xx.^3;
%  %  plot(xx,yy)
% %   ue = xx.^3;
%   
%    %error = max(error,max(abs(ue-yy)));
%    %upr = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
%    %upl = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
%    %uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
%    Z(:,i)= y;
%    
%    i=3;
% 
% y= recon4(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i+2),h(i+2),u(i+2),x(i+3),h(i+3),u(i+3),i);
% 
%   % xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
%   % yy = y(1)+y(2)*(xx-x(i))+y(3)*(xx-x(i)).^2+y(4)*(xx-x(i)).^3+y(5)*(xx-x(i)).^4;
%   % ue = exp(-(xx).^2);%xx.^3;
%    
%   % plot(xx,yy)
% %ue = xx.^3;
% 
%   % error = max(error,max(abs(ue-yy)));
%    %upr = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
%    %upl = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
%    %uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
% Z(:,i)= y;
%    
%    i = N;
%    y= recon4(x(i),h(i),u(i),x(i+1),h(i+1),u(i+1),x(i-1),h(i-1),u(i-1),x(i-2),h(i-2),u(i-2),x(i-3),h(i-3),u(i-3),i);
% 
%    %xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,M);
%    %yy = y(1)+y(2)*(xx-x(i))+y(3)*(xx-x(i)).^2+y(4)*(xx-x(i)).^3+y(5)*(xx-x(i)).^4;
%    %ue = exp(-(xx).^2);%xx.^3;
%   % plot(xx,yy)
% 
% %ue = xx.^3;
% 
% %   error = max(error,max(abs(ue-yy)))
%   %  upr = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
%   % upl = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
%   % uu(i) = u(i) + k*((upr-upl)/h(i)-f(i)); 
   
%    Z(:,i)= y;

end

