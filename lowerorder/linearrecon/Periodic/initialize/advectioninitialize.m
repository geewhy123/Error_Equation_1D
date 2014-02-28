function [u0,ue,f ] = advectioninitialize(N,x,h,tlim )
%BURGERSINITIALIZE Summary of this function goes here
%   Detailed explanation goes here
u0 = zeros(N+2,1);
ue = zeros(N+2,1);
f = zeros(N+2,1);
  for i = 2:N+1
           xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;
 

ue(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*(xr+tlim)) -cos(2*pi*(xl+tlim)));
u0(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));

f(i) = 0;
    end
f(1) = NaN;
f(N+2) = NaN;
 ue(1) = NaN;
 ue(N+2) = NaN;
 u0(1) = NaN;
 u0(N+2)= NaN;
end

