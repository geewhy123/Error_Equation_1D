function [u0,ue,f ] = poissoninitialize(N,x,h,tlim )
%BURGERSINITIALIZE Summary of this function goes here
%   Detailed explanation goes here
u0 = zeros(N+2,1);
ue = zeros(N+2,1);
f = zeros(N+2,1);
  for i = 2:N+1
         xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;   
 
       %f(i)= (1/h(i))*(1/(12*pi))*(3*cos(2*pi*xr)-cos(6*pi*xr)-3*cos(2*pi*xl)+cos(6*pi*xl));
    %f(i) = (1/h(i))*(-2*pi);%*(exp(cos(2*pi*xr))*   (2*(cos(pi*xr)^2-    4*(cos(pi*xr))^4+1)));% - (exp(cos(2*pi*xl))*(2*(cos(pi*xl))^2-4*(cos(pi*xl))^4+1));
    %f(i)=(1/h(i))*((-1/pi)*sin(2*pi*xr)+(1/(pi^3))*(sin(4*pi*xr))+(1/pi)*(sin(2*pi*xl))-(1/(pi^3))*(sin(4*pi*xl)));
   
    %%%%%f(i) = (1/h(i))*(-2*pi)*(sin(2*pi*xr)-sin(2*pi*xl));
    %f(i) =     (1/h(i))*((-4*pi^2-)/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl));


     f(i) = (1/h(i))*(pi)*(cos(pi*xr)-cos(pi*xl));%(1/h(i))*(-4*pi^2)*( (exp(1)^3*sin(2*pi*xr)+1)/(sin(2*pi*xr)+exp(1)^3)^2 - (exp(1)^3*sin(2*pi*xl)+1)/(sin(2*pi*xl)+exp(1)^3)^2);
   %f(i) = -(1/h(i))*(2*pi)*(sin(2*pi*xr)-sin(2*pi*xl));
    
  
    
    %u(i) = (1/h(i))*((x(i)+h(i)/2)^4-((x(i)-h(i)/2)^4))/4;%exp(-(x(i)-0.5)^2);
    %u(i) = (1/h(i))*((x(i)+h(i)/2)^5-((x(i)-h(i)/2)^5))/5;%exp(-(x(i)-0.5)^2);
    %u(i) = (1/h(i))*((pi^(1/2)*erf(x(i)+h(i)/2 ))/2-(pi^(1/2)*erf(x(i)-h(i)/2 ))/2);
    %u(i) = (1/h(i))*(-1/(2*pi))*(exp(cos(2*pi*(x(i)+h(i)/2)))-exp(cos(2*pi*(x(i)-h(i)/2))));    
    %u(i) = (1/h(i))*(-pi)*(cos(2*pi*(x(i)+h(i)/2))-3*cos(6*pi*(x(i)+h(i)/2))+cos(2*pi*(x(i)-h(i)/2))+3*cos(6*pi*(x(i)-h(i)/2)));
 
   
    %u(i) = (1/h(i))*(4*pi*sin(2*pi*xr)-(16/pi)*sin(4*pi*xr)-4*pi*sin(2*pi*xl)+(16/pi)*sin(4*pi*xl));
    
%%%%%    ue(i) = (1/h(i))*(1/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl));


 %this%% ue(i) = (1/h(i))*(log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl)));
 
 ue(i)= (1/h(i))*(1/pi)*(-cos(pi*xr)+cos(pi*xl));%(1/h(i))*((-1/(2*pi))*(100*exp(-4*pi^2*tlim))*(cos(2*pi*xr)-cos(2*pi*xl))+  (log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl))));%(1/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl)));



 %initial
 u0(i) = 0;%(1/h(i))*((-1/(2*pi))*(100*(cos(2*pi*xr)-cos(2*pi*xl))) +(log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl))));%(1/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl)));
 
%u0(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
    
    end
f(1) = NaN;
f(N+2) = NaN;
 ue(1) = NaN;
 ue(N+2) = NaN;
 u0(1) = NaN;
 u0(N+2)= NaN;
end

