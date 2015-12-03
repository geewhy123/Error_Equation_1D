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
elseif(strcmp(obj.physics,'BurgersVisc')==1)
    burgersviscinitialize(obj);
elseif(strcmp(obj.physics,'EulerQ')==1)
    eulerinitialize(obj);
elseif(strcmp(obj.physics,'LinearSystem')==1)
    linearsysteminitialize(obj);
else
    fprintf('No exact solution provided')
%     error('1')
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
une = zeros(N+2,1);
ne = zeros(N+2,1);
f = zeros(N+2,1);
z = rand(N+2,1);
% z = [0; randn(N,1); 0];
rr  = randn(N+2,1);
rr = rr - sum(rr(2:N+1)/N);
s= 0;
  for i = 2:N+1
         xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;   
 
       %f(i)= (1/h(i))*(1/(12*pi))*(3*cos(2*pi*xr)-cos(6*pi*xr)-3*cos(2*pi*xl)+cos(6*pi*xl));
    %f(i) = (1/h(i))*(-2*pi);%*(exp(cos(2*pi*xr))*   (2*(cos(pi*xr)^2-    4*(cos(pi*xr))^4+1)));% - (exp(cos(2*pi*xl))*(2*(cos(pi*xl))^2-4*(cos(pi*xl))^4+1));
    %f(i)=(1/h(i))*((-1/pi)*sin(2*pi*xr)+(1/(pi^3))*(sin(4*pi*xr))+(1/pi)*(sin(2*pi*xl))-(1/(pi^3))*(sin(4*pi*xl)));
   
    %%%%%f(i) = (1/h(i))*(-2*pi)*(sin(2*pi*xr)-sin(2*pi*xl));
    %f(i) =     (1/h(i))*((-4*pi^2-)/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl));

if(obj.bcLeftType == 'P' && obj.bcRightType == 'P' || strcmp(obj.goal,'SS')~=1)
   
%      f(i) = (1/h(i))*(-4*pi^2)*( (exp(1)^3*sin(2*pi*xr)+1)/(sin(2*pi*xr)+exp(1)^3)^2 - (exp(1)^3*sin(2*pi*xl)+1)/(sin(2*pi*xl)+exp(1)^3)^2);
     f(i) = (1/h(i))*(2*pi)*(cos(2*pi*xr)-cos(2*pi*xl));
%      z(i) = randn();
%      f(i) = f(i) + (1/h(i))*z(i)*1e-1*(xr-xl);
elseif(obj.bcLeftType == 'D' && obj.bcRightType == 'D')
f(i) = (1/h(i))*pi*(cos(pi*xr)-cos(pi*xl));
% f(i) = (1/h(i))*2*pi*(cos(2*pi*xr)-cos(2*pi*xl));
% f(i) = (1/h(i))*(2*xr-3*xr^2-2*xl+3*xl^2);



% f(i) = (1/h(i))*(-0.5)*(sech(xr/2)^2-sech(xl/2)^2);
% f(i) = (1/h(i))*pi*(cos(pi*(xr+0.2))-cos(pi*(xl+0.2)));
f(i) = (1/h(i))*2*pi*(cos(2*pi*xr)-cos(2*pi*xl));
m = [9 13];%rand();%randn();
% m = [2:20];
m = 5;
for j = 1:length(m)*0
f(i) = f(i) + rr(i)*1e-2*(1/h(i))*2*m(j)*pi*(cos(2*pi*m(j)*xr)-cos(2*pi*m(j)*xl));
end

% % % f(i) = f(i) + 1e-2*(1/h(i))*2*3*pi*(cos(30*pi*(xr+0.312))-cos(pi*30*(xl+0.312))) - 1e-2*(1/h(i))*2*5*pi*(cos(2*pi*5*(xr-0.5))-cos(2*pi*5*(xl-0.5)));

z(i) = randn();

f(i) = f(i) + (1/h(i))*z(i)*1e-2*(xr-xl)*0;

% if(i == N+1)
%    f(i) = f(i) - z(i)*1e1 ;
%    z(N+1) = -sum(z(2:N));
%    f(i) = f(i) + z(i)*1e1;
% end

s = s+f(i);

% if(i==N+1)
% f = f-s/N;    
% % z = z-sum/N;
% end



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
 if(obj.bcLeftType == 'P' && obj.bcRightType == 'P'|| strcmp(obj.goal,'SS')~=1)

ue(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
if(strcmp(obj.goal,'TimeAccurate')==1)
 ue(i)= (1/h(i))*((-1/(2*pi))*(100*exp(-4*pi^2*tlim))*(cos(2*pi*xr)-cos(2*pi*xl)));%+  (log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl))));%(1/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl)));    
%   ue(i)= (1/h(i))*((-1/(pi))*(100*exp(-pi^2*tlim))*(cos(pi*xr)-cos(pi*xl)));
 f(i) = 0;
end


 elseif(obj.bcLeftType == 'D' && obj.bcRightType == 'D')
 ue(i) = (1/h(i))*((1/pi)*(-cos(pi*xr)+cos(pi*xl)));
%  ue(i) = (1/h(i))*((1/(2*pi))*(-cos(2*pi*xr)+cos(2*pi*xl)));
 
%  ue(i) = (1/h(i))*(xr^3/3-xr^4/4-xl^3/3+xl^4/4);
% ue(i) = (1/h(i))*(0.925*(xr-xl)+(0.025/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl)));
%  ue(i) = (1/h(i))*(xr-2*log(cosh((xr)/2))-xl+2*log(cosh((xl)/2)));
% ue(i) = (1/h(i))*((1/pi))*(-cos(pi*(xr+0.2))+cos(pi*(xl+0.2)));
ue(i) = (1/h(i))*((1/(2*pi))*(-cos(2*pi*xr)+cos(2*pi*xl)));


ne(i) =  (1/h(i))*z(i)*1e1*(xr^3/6-xl^3/6);

if(i==N+1)
   une = ue+ ne-s/N;  
end
 else
    assert(0) 
 end
% ue(i) = (1/h(i))*((xr-0.5)^5/5-(xl-0.5)^5/5);

 %initial
% % %  u0(i) = (1/h(i))*((-1/(2*pi))*(100*(cos(2*pi*xr)-cos(2*pi*xl))) +(log(exp(1)^3+sin(2*pi*xr))-log(exp(1)^3+sin(2*pi*xl))));%(1/(2*pi))*(sin(2*pi*xr)-sin(2*pi*xl)));
 u0(i)=0;
 
 if(strcmp(obj.goal,'TimeAccurate')==1)
 u0(i)=  (1/h(i))*(-100/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
%   u0(i)=  (1/h(i))*(-100/(pi))*(cos(pi*xr)-cos(pi*xl));
end
% u0(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
    
  end
%     plot(x,u0,x,ue)
%     error('1')
f(1) = NaN;
f(N+2) = NaN;
 ue(1) = NaN;
 ue(N+2) = NaN;
 u0(1) = NaN;
 u0(N+2)= NaN;
 

 obj.exactSolution = ue;
 obj.initialSolution = u0;
 obj.source = f;
 obj.exactNoisySolution = une*0;
 

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
if(strcmp(obj.goal,'SS')==1)%|| obj.endTime > 8)
   ue(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
   u0(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
   f(i) = (1/h(i))*(sin(2*pi*xr)-sin(2*pi*xl)); 
end
%  ue(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
%    u0(i) = (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
% f(i) = (1/h(i))*(sin(2*pi*xr)-sin(2*pi*xl)); 

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
     

     
%     MMS
%     ue(i) = (1/h(i))*((1/pi)*(-cos(pi*xr)+cos(pi*xl)));
%     u0 = ue;
%     f(i) = -(1/h(i))*( ( - sin(pi*xr) - pi*cos(pi*xr) - cos(pi*xr)^2/2) -(- sin(pi*xl) - pi*cos(pi*xl) - cos(pi*xl)^2/2));
else
    assert(0)
end

% f(i) = 0;
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

function  burgersviscinitialize(obj)
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
    ue(i) = (1/h(i))*((1/(2*pi))*(-cos(2*pi*xr)+cos(2*pi*xl)));
     u0(i) = ue(i);
%     f(i) = -(2/h(i))*((2*pi+sin(pi*xr)^2)*(sin(pi*xr)^2-1) - (2*pi+sin(pi*xl)^2)*(sin(pi*xl)^2-1));
% f(i) = -(1/h(i))*( (2*(2*pi + sin(pi*xr)^2)*(sin(pi*xr)^2 - 1)) - (2*(2*pi + sin(pi*xl)^2)*(sin(pi*xl)^2 - 1)) )
f(i) = (-2/h(i))*(  ((sin(pi*xr)^2 - 1)*(2*pi - sin(pi*xr)^2)) - ((sin(pi*xl)^2 - 1)*(2*pi - sin(pi*xl)^2)));
 elseif(obj.bcLeftType == 'D' && obj.bcRightType == 'D')
    
    ue(i) = (1/h(i))*-2*(log(cosh(xr))-log(cosh(xl)));
     u0(i) = (1/h(i))*-2*(log(cosh(xr))-log(cosh(xl)));
     

     if(strcmp(obj.goal,'TimeAccurate')==1)
        b = 10;
        u0(i) = u0(i)+ (1/h(i))*((xr^3/3-xr^2/2)-(xl^3/3-xl^2/2));
        ue(i) = ue(i) + (1/h(i))*exp(-b*obj.endTime)*((xr^3/3-xr^2/2)-(xl^3/3-xl^2/2));
%         f(i) = (1/h(i))*(-b*exp(-b*obj.endTime)*((xr^3/3-xr^2/2)-(xl^3/3-xl^2/2)) + exp(-2*b*obj.endTime)*((xr^4/2-xr^3+xr^2/2)-(xl^4/2-xl^3+xl^2/2)) );
     end

     
%     MMS
%     ue(i) = (1/h(i))*((1/pi)*(-cos(pi*xr)+cos(pi*xl)));
%     u0 = ue;
%     f(i) = -(1/h(i))*( ( - sin(pi*xr) - pi*cos(pi*xr) - cos(pi*xr)^2/2) -(- sin(pi*xl) - pi*cos(pi*xl) - cos(pi*xl)^2/2));
else
    assert(0)
end

% f(i) = 0;
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
a = 0.0;
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
%    c1 = 0.3478548451;
% c2 = 0.6521451549;
% c3 = 0.6521451549;
% c4 = 0.3478548451;
% x1= 0.8611363116;
% x2 = 0.339981436;
% x3 = -0.339981436;
% x4= -0.8611363116;
% 
%  xx1 = ((xr-xl)/2)*x1+(xr+xl)/2;
%  xx2 = ((xr-xl)/2)*x2+(xr+xl)/2;
%  xx3 = ((xr-xl)/2)*x3+(xr+xl)/2;
%  xx4 = ((xr-xl)/2)*x4+(xr+xl)/2;


%  cc = [0.3478548451;
%  0.6521451549;
%  0.6521451549;
%  0.3478548451];
% yy = [ 0.8611363116;
%  0.339981436;
%  -0.339981436;
%  -0.8611363116];
cc = [0.1713244924
    0.3607615730
    0.4679139346
    0.4679139346
    0.3607615730
    0.1713244924
    ];
yy = [0.9324695142
    0.6612093865
    0.2386191861
    -0.2386191861
    -0.6612093865
    -0.9324695142
    ];

xx = ((xr-xl)/2)*yy+(xr+xl)/2;



ts = 4;%/(2*pi);
A = (1/(2*pi*ts));

for j = 1:length(xx);
% F = @(s) s+tlim*(sin(2*pi*s)*A)-xx1;
% xxx1=fzero(F,0);
% F = @(s) s+tlim*(sin(2*pi*s)*A)-xx2;
% xxx2=fzero(F,0);
% F = @(s) s+tlim*(sin(2*pi*s)*A)-xx3;
% xxx3=fzero(F,0);
% F = @(s) s+tlim*(sin(2*pi*s)*A)-xx4;
% xxx4=fzero(F,0);
F = @(s) s+tlim*(sin(2*pi*s)*A)-xx(j);
xxx(j)=fzero(F,0);
end

% ue(i) = a+ (1/(2*pi))* (1/h(i))*((xr-xl)/2)*(c1*sin(2*pi*xxx1)+c2*sin(2*pi*xxx2)+c3*sin(2*pi*xxx3)+c4*sin(2*pi*xxx4));
ue(i) = A* (1/h(i))*((xr-xl)/2)*dot(cc,sin(2*pi*xxx));%(c1*sin(2*pi*xxx1)+c2*sin(2*pi*xxx2)+c3*sin(2*pi*xxx3)+c4*sin(2*pi*xxx4));
% ue(i) = (1/(2*pi))* (1/h(i))*((xr-xl)/2)*(c1*(sin(2*pi*xxx1)+4*pi*h(i))+c2*(sin(2*pi*xxx2)+4*pi*h(i))+c3*(sin(2*pi*xxx3)+4*pi*h(i))+c4*(sin(2*pi*xxx4)+4*pi*h(i)));

% u0(i) = a+(1/(2*pi))* (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
u0(i) = -A* (1/h(i))*(-1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));


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

ur = 0;
ul = 0;

  if(tlim < 1e-12)
     ur  = cos(2*pi*xr);
     ul = cos(2*pi*xl);
     ue(i) = (1/h(i))*(1/(2*pi))*(ur-ul);
  else
  for n = 1:100
     ur = ur +besselj(n,n*tlim)*-cos(n*2*pi*xr)/(n*tlim*n);
     ul = ul +besselj(n,n*tlim)*-cos(n*2*pi*xl)/(n*tlim*n);
  end
  
  ue(i) = (-2/(2*pi))*(1/h(i))*(1/(2*pi))*(ur-ul);
  end
  
u0(i) = (1/(2*pi))*(1/h(i))*(1/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));



  end

% ue
% 
% %   
%     x = linspace(0,1,N+2)';
% 
%   u = 0*x;
% 
%   if(tlim < 1e-12)
%      u  = -sin(x);
%   else
%   for n = 1:100
%      u = u+besselj(n,n*tlim)*sin(n*2*pi*x)/(n*tlim);
%   end
%   u = -2*u/(2*pi);
%   end
  
%   

%   ue =u ;
%   u0 = -sin(2*pi*x)/(2*pi);
% u
% error('1')
  
  
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





function eulerinitialize(obj)
x = obj.cellCentroids;
N = obj.nCells;
h = obj.cellWidths;
U = NaN*ones(N+2,3);


 for i = 2:N+1
           xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;

% U(i,1) = (1/h(i))*(0.5*(xr-xl)+(1/pi)*(sin(pi*xr)-sin(pi*xl)));
% U(i,2) = (1/h(i))*(1*(xr-xl)-(1/pi)*(cos(pi*xr)-cos(pi*xl)));
% U(i,3) = (1/h(i))*(0.4*(xr-xl)+(1/pi)*(sin(pi*xr)-sin(pi*xl)));

UL = obj.bcLeftVal;
UR = obj.bcRightVal;

% a = asin((UL(1)-1)/.5);
% b = asin((UR(1)-1)/.5)-a;
% U(i,1) = (1/h(i))*(xr-(0.5/a)*cos(a*xr+b)-xl+(0.5/a)*cos(a*xl+b));
% a = asin((UL(2)-1)/.5);
% b = asin((UR(2)-1)/.5)-a;
% U(i,2) = (1/h(i))*(xr-(0.5/a)*cos(a*xr+b)-xl+(0.5/a)*cos(a*xl+b));
% a = asin((UL(3)-1)/.5);
% b = asin((UR(3)-1)/.5)-a;
% U(i,3) = (1/h(i))*(xr-(0.5/a)*cos(a*xr+b)-xl+(0.5/a)*cos(a*xl+b));

% % % a = UR(1)-UL(1)-1;
% % % U(i,1) = (1/h(i))*(UL(1)*xr+(a/2)*xr^2+xr^3/3-UL(1)*xl-(a/2)*xl^2-xl^3/3);
% % % a = UR(2)-UL(2)-1;
% % % U(i,2) = (1/h(i))*(UL(2)*xr+(a/2)*xr^2+xr^3/3-UL(2)*xl-(a/2)*xl^2-xl^3/3);
% % % a = UR(3)-UL(3)-1;
% % % U(i,3) = (1/h(i))*(UL(3)*xr+(a/2)*xr^2+xr^3/3-UL(3)*xl-(a/2)*xl^2-xl^3/3);
%     gam = 1.4;
% %     Pa = obj.bcLeftVal(3);
% Pa = 0.95
%     T0=1;
%     P0=1;
%     Ta = T0*(Pa/P0)^((gam-1)/gam);
%     rhoa = Pa/Ta;%*gam
%     ua = sqrt((2/(gam-1))*(T0/Ta-1)) *sqrt(gam*Pa/rhoa);%???
%     fprintf('check quant definitions, and consistent')


U(i,1) = (1/h(i))*0.964024984265400*(xr-xl)+(1/h(i))*0.0*(xr-xl) + (1/h(i))*(0.01/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
U(i,2) = (1/h(i))*0.319121501427128*(xr-xl)+(1/h(i))*0.0*(xr-xl) + (1/h(i))*(0.01/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));
U(i,3) = (1/h(i))*0.95*(xr-xl)             +(1/h(i))*0.0*(xr-xl) + (1/h(i))*(0.01/(2*pi))*(cos(2*pi*xr)-cos(2*pi*xl));




%%%

% pleft = 1.0;  pright = 0.1; rholeft = 1.0;  rhoright = 0.125;
% uleft = 0;  uright = 0; tend = 0.17; lambda = 0.35;
%  if x(i) < 0.5, U(i,3) = (1/h(i))*pleft*(xr-xl); U(i,1) = (1/h(i))*rholeft*(xr-xl);  U(i,2) = (1/h(i))*uleft*(xr-xl); 
%  else,	      U(i,3) = (1/h(i))*pright*(xr-xl); U(i,1) = (1/h(i))*rhoright*(xr-xl); U(i,2) = (1/h(i))*uright*(xr-xl);
%  end 

  
%%%


f(i,1:3) = 0;
    end

obj.initialSolution = U;


end

function  linearsysteminitialize(obj)
%BURGERSINITIALIZE Summary of this function goes here
%   Detailed explanation goes here
tlim = obj.endTime;
x = obj.cellCentroids;
N = obj.nCells;
h = obj.cellWidths;

u0 = zeros(N+2,2);
ue = zeros(N+2,2);
f = zeros(N+2,1);
  for i = 2:N+1
           xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;
 
    
    if(obj.bcLeftType == 'D' && obj.bcRightType == 'D')
        ue(i,1) = (1/h(i))*( (-(10^(1/2)*pi^(1/2)*(erf(10^(1/2)*(10*tlim - 10*xr + 5)) - erf(10^(1/2)*(10*tlim + 10*xr - 5))))/400) - ...
                           (-(10^(1/2)*pi^(1/2)*(erf(10^(1/2)*(10*tlim - 10*xl + 5)) - erf(10^(1/2)*(10*tlim + 10*xl - 5))))/400) );
                       
                       ue(i,2)=(1/h(i))*( (  ( exp(-1000*(xr-tlim-0.5).^2) + exp(-1000*(xr+tlim-0.5).^2))/2) - ...
                                          (  ( exp(-1000*(xl-tlim-0.5).^2) + exp(-1000*(xl+tlim-0.5).^2))/2) );
                       
                       
                       u0(i,1) = (1/h(i))*( (-(10^(1/2)*pi^(1/2)*(erf(10^(1/2)*( - 10*xr + 5)) - erf(10^(1/2)*( + 10*xr - 5))))/400) - ...
                           (-(10^(1/2)*pi^(1/2)*(erf(10^(1/2)*( - 10*xl + 5)) - erf(10^(1/2)*( + 10*xl - 5))))/400) );
        
                        u0(i,2)=(1/h(i))*( (  ( exp(-1000*(xr-0.5).^2) + exp(-1000*(xr-0.5).^2))/2) - ...
                                          (  ( exp(-1000*(xl-0.5).^2) + exp(-1000*(xl-0.5).^2))/2) );
                       
% if((xr+xl)/2 > 0.375 && (xr+xl)/2 < 0.625)
% u0(i,1) = (1/h(i))*((1/(8*pi))*(sin(8*pi*(xr-0.5))-sin(8*pi*(xl-0.5))) + xr-xl);
% u0(i,2) = (1/h(i))*(cos(8*pi*(xr-0.5))-cos(8*pi*(xl-0.5)) + xr-xl);
% else 
%     u0(i,1) = 0;
%     u0(i,2) = 0;
% end
% 
% if((xr+xl)/2 > 0.375-tlim && (xr+xl)/2 < 0.625+tlim)
% ue(i,1) = ((1/h(i))*((1/(8*pi))*(sin(8*pi*(xr+tlim-0.5))-sin(8*pi*(xl+tlim-0.5))) + xr+tlim-(xl+tlim)) + ...
%            (1/h(i))*((1/(8*pi))*(sin(8*pi*(xr-tlim-0.5))-sin(8*pi*(xl-tlim-0.5))) + xr-tlim-(xl-tlim)))/2;
% end
% ue

    elseif(obj.bcLeftType == 'P' && obj.bcRightType == 'P')
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

