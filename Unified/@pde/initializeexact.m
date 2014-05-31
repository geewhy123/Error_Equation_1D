function  initializeexact(obj )
%INITIALIZEEXACT Summary of this function goes here
%   Detailed explanation goes here

if(strcmp(obj.physics,'Burgers')==1)
   burgersinitialize(obj);

elseif(strcmp(obj.physics,'Advection')==1)
     advectioninitialize(obj);
elseif(strcmp(obj.physics,'Poisson')==1)
poissoninitialize(obj);
elseif(strcmp(obj.physics,'BurgersMod')==1)
    burgersmodinitialize(obj);
else
    error('1')
end


end

function poissoninitialize(obj)
%    function poissoninitialize(problem )
%BURGERSINITIALIZE Summary of this function goes here
%   Detailed explanation goes here
tlim = obj.endTime;
x = obj.cellCentroids;
N = obj.nCells;
h = obj.cellWidths;

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

if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')
     f(i) = (1/h(i))*(-4*pi^2)*( (exp(1)^3*sin(2*pi*xr)+1)/(sin(2*pi*xr)+exp(1)^3)^2 - (exp(1)^3*sin(2*pi*xl)+1)/(sin(2*pi*xl)+exp(1)^3)^2);
%      f(i) = (1/h(i))*(2*pi)*(cos(2*pi*xr)-cos(2*pi*xl));
elseif(obj.bcLeftType == 'D' && obj.bcRightType == 'D')
f(i) = (1/h(i))*pi*(cos(pi*xr)-cos(pi*xl));
% f(i) = (1/h(i))*(2*xr-3*xr^2-2*xl+3*xl^2);
else
   assert(0) 
end
%    f(i) = -(1/h(i))*(2*pi)*(sin(2*pi*xr)-sin(2*pi*xl));
    
  
    
    %u(i) = (1/h(i))*((x(i)+h(i)/2)^4-((x(i)-h(i)/2)^4))/4;%exp(-(x(i)-0.5)^2);
    %u(i) = (1/h(i))*((x(i)+h(i)/2)^5-((x(i)-h(i)/2)^5))/5;%exp(-(x(i)-0.5)^2);
    %u(i) = (1/h(i))*((pi^(1/2)*erf(x(i)+h(i)/2 ))/2-(pi^(1/2)*erf(x(i)-h(i)/2 ))/2);
    %u(i) = (1/h(i))*(-1/(2*pi))*(exp(cos(2*pi*(x(i)+h(i)/2)))-exp(cos(2*pi*(x(i)-h(i)/2))));    
    %u(i) = (1/h(i))*(-pi)*(cos(2*pi*(x(i)+h(i)/2))-3*cos(6*pi*(x(i)+h(i)/2))+cos(2*pi*(x(i)-h(i)/2))+3*cos(6*pi*(x(i)-h(i)/2)));
 
   
    %u(i) = (1/h(i))*(4*pi*sin(2*pi*xr)-(16/pi)*sin(4*pi*xr)-4*pi*sin(2*pi*xl)+(16/pi)*sin(4*pi*xl));
    
%%%%%    ue(i) = (1/h(i))*(1/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl));


 %this%% ue(i) = (1/h(i))*(log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl)));
 if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')
 ue(i)= (1/h(i))*((-1/(2*pi))*(100*exp(-4*pi^2*tlim))*(cos(2*pi*xr)-cos(2*pi*xl))+  (log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl))));%(1/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl)));
%  ue(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
 elseif(obj.bcLeftType == 'D' && obj.bcRightType == 'D')
 ue(i) = (1/h(i))*((1/pi)*(-cos(pi*xr)+cos(pi*xl)));
%  ue(i) = (1/h(i))*(xr^3/3-xr^4/4-xl^3/3+xl^4/4);
 else
    assert(0) 
 end
% ue(i) = (1/h(i))*((xr-0.5)^5/5-(xl-0.5)^5/5);

 %initial
% % %  u0(i) = (1/h(i))*((-1/(2*pi))*(100*(cos(2*pi*xr)-cos(2*pi*xl))) +(log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl))));%(1/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl)));
 u0(i)=0;
% u0(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
    
    end
f(1) = NaN;
f(N+2) = NaN;
 ue(1) = NaN;
 ue(N+2) = NaN;
 u0(1) = NaN;
 u0(N+2)= NaN;
 
 
 obj.exactSolution = ue;
 obj.initialSolution = u0;
 obj.source = f;
 

end


function  advectioninitialize(obj)
%BURGERSINITIALIZE Summary of this function goes here
%   Detailed explanation goes here
tlim = obj.endTime;
x = obj.cellCentroids;
N = obj.nCells;
h = obj.cellWidths;

u0 = zeros(N+2,1);
ue = zeros(N+2,1);
f = zeros(N+2,1);
  for i = 2:N+1
           xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;
 
if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')
ue(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*(xr+tlim)) -cos(2*pi*(xl+tlim)));
u0(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
elseif(obj.bcLeftType == 'F' && obj.bcRightType == 'D')
    ue(i) = (1/h(i))*(xr-xl);%(1/h(i))*(1/pi)*(-cos(pi*xr)+cos(pi*xl));%(1/h(i))*(xr-xl);
    u0(i) = 0;%(1/h(i))*((1/pi)*(-cos(pi*xr)+cos(pi*xl))+xr-xl);
else
    assert(0)
end

f(i) = 0;
    end
f(1) = NaN;
f(N+2) = NaN;
 ue(1) = NaN;
 ue(N+2) = NaN;
 u0(1) = NaN;
 u0(N+2)= NaN;
 obj.exactSolution = ue;
 obj.initialSolution = u0;
 obj.source = f;
end

function  burgersmodinitialize(obj)
%BURGERSINITIALIZE Summary of this function goes here
%   Detailed explanation goes here
tlim = obj.endTime;
x = obj.cellCentroids;
N = obj.nCells;
h = obj.cellWidths;

u0 = zeros(N+2,1);
ue = zeros(N+2,1);
f = zeros(N+2,1);
  for i = 2:N+1
           xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;
 
if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')
ue(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*(xr+tlim)) -cos(2*pi*(xl+tlim)));
u0(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
elseif(obj.bcLeftType == 'F' && obj.bcRightType == 'D')
    ue(i) = (1/h(i))*(xr-2*log(cosh((xr)/2))-xl+2*log(cosh((xl)/2)));
    u0(i) = (1/h(i))*(xr-2*log(cosh((xr)/2))-xl+2*log(cosh((xl)/2)));
elseif(obj.bcLeftType == 'D' && obj.bcRightType == 'D')
    ue(i) = (1/h(i))*(xr-2*log(cosh((xr)/2))-xl+2*log(cosh((xl)/2)));
     u0(i) = (1/h(i))*(xr-2*log(cosh((xr)/2))-xl+2*log(cosh((xl)/2)));%(1/h(i))*(xr-tanh(1/2)*xr^2/2-xl+tanh(1/2)*xl^2);%
else
    assert(0)
end

f(i) = 0;
    end
f(1) = NaN;
f(N+2) = NaN;
 ue(1) = NaN;
 ue(N+2) = NaN;
 u0(1) = NaN;
 u0(N+2)= NaN;
 obj.exactSolution = ue;
 obj.initialSolution = u0;
 obj.source = f;
end


function  burgersinitialize(obj )
%BURGERSINITIALIZE Summary of this function goes here
%   Detailed explanation goes here
a = 0.5;
% for k = 1:30
% tlim = tlim+0.01;
    

tlim = obj.endTime;
x = obj.cellCentroids;
N = obj.nCells;
h = obj.cellWidths;
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
 
 
  obj.exactSolution = ue;
 obj.initialSolution = u0;
 obj.source = f;
end

% tlim

% end



