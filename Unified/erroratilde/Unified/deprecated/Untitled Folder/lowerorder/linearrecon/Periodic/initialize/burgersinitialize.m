function [u0,ue,f ] = burgersinitialize(N,x,h,tlim )
%BURGERSINITIALIZE Summary of this function goes here
%   Detailed explanation goes here
a = 0.5;
% for k = 1:30
% tlim = tlim+0.01;
    
u0 = NaN*zeros(N+2,1);
ue = NaN*zeros(N+2,1);
ue2 = NaN*zeros(N+2,1);

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


F = @(s) s+tlim*(sin(2*pi*s)/(2*pi)+a)-xx1;
xxx1=fzero(F,0);
F = @(s) s+tlim*(sin(2*pi*s)/(2*pi)+a)-xx2;
xxx2=fzero(F,0);
F = @(s) s+tlim*(sin(2*pi*s)/(2*pi)+a)-xx3;
xxx3=fzero(F,0);
F = @(s) s+tlim*(sin(2*pi*s)/(2*pi)+a)-xx4;
xxx4=fzero(F,0);


ue(i) = a+ (1/(2*pi))* (1/h(i))*((xr-xl)/2)*(c1*sin(2*pi*xxx1)+c2*sin(2*pi*xxx2)+c3*sin(2*pi*xxx3)+c4*sin(2*pi*xxx4));
% ue(i) = (1/(2*pi))* (1/h(i))*((xr-xl)/2)*(c1*(sin(2*pi*xxx1)+4*pi*h(i))+c2*(sin(2*pi*xxx2)+4*pi*h(i))+c3*(sin(2*pi*xxx3)+4*pi*h(i))+c4*(sin(2*pi*xxx4)+4*pi*h(i)));

u0(i) = a+(1/(2*pi))* (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));


f(i) = 0;


%new


% ucur = (1/(2*pi))*sin(2*pi*(xx1))+a;
% unext = 123;
% while(abs(unext-ucur)>1e-15)
%    ucur = unext;
%    unext = (1/(2*pi))*sin(2*pi*(xx1-ucur*tlim))+a   ;
% end
% u1 = unext;


% ucur = (1/(2*pi))*sin(2*pi*(xx2))+a;
% unext = 123;
% while(abs(unext-ucur)>1e-15)
%    ucur = unext;
%    unext = (1/(2*pi))*sin(2*pi*(xx2-ucur*tlim))+a;   
% end
% u2 = unext;
% 
% ucur = (1/(2*pi))*sin(2*pi*(xx3))+a;
% unext = 123;
% while(abs(unext-ucur)>1e-15)
%    ucur = unext;
%    unext = (1/(2*pi))*sin(2*pi*(xx3-ucur*tlim))+a;   
% end
% u3 = unext;
% 
% ucur = (1/(2*pi))*sin(2*pi*(xx4))+a;
% unext = 123;
% while(abs(unext-ucur)>1e-15)
%    ucur = unext;
%    unext = (1/(2*pi))*sin(2*pi*(xx4-ucur*tlim))+a;   
% end
% u4 = unext;
% 
% ue2(i) = (1/h(i))*((xr-xl)/2)*(c1*u1+c2*u2+c3*u3+c4*u4);


%new

  end
  
%   
% plot(x,u0,'*-',x,ue,'+-')
% error('1')

% max(abs(ue-ue2))
%   plot(x,u0,'*-',x,ue,'+-',x,ue2,'x-')
%   error('1')
% D(k) = getframe;


f(1) = NaN;
f(N+2) = NaN;
 ue(1) = NaN;
 ue(N+2) = NaN;
 u0(1) = NaN;
 u0(N+2)= NaN;
end

% tlim

% end

