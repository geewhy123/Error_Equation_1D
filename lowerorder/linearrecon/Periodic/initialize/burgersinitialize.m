function [u0,ue,f ] = burgersinitialize(N,x,h,tlim )
%BURGERSINITIALIZE Summary of this function goes here
%   Detailed explanation goes here
u0 = zeros(N+2,1);
ue = zeros(N+2,1);
f = zeros(N+2,1);
  for i = 2:N+1
           xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;
   c1 = 0.3478548451;
c2 = 0.6521451549;
c3 = 0.6521451549;
c4 = 0.3478548451;
x1= 0.8611363116;
x2 = 0.339981436;
x3 = -0.339981436;
x4= -0.8611363116;

 xx1 = ((xr-xl)/2)*x1+(xr+xl)/2;
 xx2 = ((xr-xl)/2)*x2+(xr+xl)/2;
 xx3 = ((xr-xl)/2)*x3+(xr+xl)/2;
 xx4 = ((xr-xl)/2)*x4+(xr+xl)/2;


F = @(s) s+tlim*sin(2*pi*s)/(2*pi)-xx1;
xxx1=fzero(F,0);
F = @(s) s+tlim*sin(2*pi*s)/(2*pi)-xx2;
xxx2=fzero(F,0);
F = @(s) s+tlim*sin(2*pi*s)/(2*pi)-xx3;
xxx3=fzero(F,0);
F = @(s) s+tlim*sin(2*pi*s)/(2*pi)-xx4;
xxx4=fzero(F,0);


ue(i) = (1/(2*pi))* (1/h(i))*((xr-xl)/2)*(c1*sin(2*pi*xxx1)+c2*sin(2*pi*xxx2)+c3*sin(2*pi*xxx3)+c4*sin(2*pi*xxx4));

u0(i) = (1/(2*pi))* (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));

f(i) = 0;
    end
f(1) = NaN;
f(N+2) = NaN;
 ue(1) = NaN;
 ue(N+2) = NaN;
 u0(1) = NaN;
 u0(N+2)= NaN;
end

