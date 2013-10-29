clear all
close all
N =10;
h = 1/N;
A = zeros(N+2,N+2);

x = zeros(N+2,1);
for i = 1:N+2
   x(i) = (i-1.5)*h; 
end

for i = 1:N+2
   for j = 1:N+2
      if(i==j)
         A(i,i) = -2/h^2; 
      end
      if(i+1==j || i-1==j)
         A(i,j) = 1/h^2; 
      end
   end
end

A(1,1) = 1;
A(1,2) = 1;
A(N+2,N+1)=1;
A(N+2,N+2)=1;

%A(1,N) = 1;
%A(N,1) = 1;

f = zeros(N+2,1);
for i = 1:N+2
    %f(i) = exp(-x(i)^2)*(-2+4*x(i)^2); 
    f(i) = (1/h)*(-2*((i-1)*h)*exp(-((i-1)*h)^2)+2*((i-2)*h)*exp(-((i-2)*h)^2));
end

f(1)=2*1;
f(N+2)=2*exp(-1);
u = A\f;



E = zeros(N+2,1);
s = zeros(N+2,1);
for i = 1:N+2
    if ((i > 2) && (i < N+1))
    s(i) =  (-u(i-2)+16*u(i-1)-30*u(i)+16*u(i+1)-u(i+2))/(12*h^2);
    end
    %if ((i > 1) && (i < N+2))
    %s(i) =  (u(i-1)-2*u(i)+u(i+1))/(h^2);
    %end
    
   % s(i) = (-2+4*((i-1.5)*h)^2)*(exp(-((i-1.5)*h)^2));
    s(i) = s(i) - (1/h)*(-2*((i-1)*h)*exp(-((i-1)*h)^2)+2*((i-2)*h)*exp(-((i-2)*h)^2));
   
end
s(2)=s(2)+(10*u(1)-15*u(2)-4*u(3)+14*u(4)-6*u(5)+u(6))/(12*h^2);

s(N+1)=s(N+1)+(10*u(N+2)-15*u(N+1)-4*u(N)+14*u(N-1)-6*u(N-2)+u(N-3))/(12*h^2);



s = -s;
%s(1)=2*(1-(u(1)+u(2))/2);
%s(N+2)=2*(exp(-1)-(u(N+1)+u(N+2))/2);
%s = -s;
%E=A\s;





uea = zeros(N+2,1);
for i = 1:N+2
uea(i) = (1/h)*(sqrt(pi)*erf((i-1)*h)/2-sqrt(pi)*erf((i-2)*h)/2);
end

%s(1) = 2*(uea(1)-u(1));
%s(N+2)=2*(uea(N+2)-u(N+2));
s(1) = 2*(1-(25*u(2)-23*u(3)+13*u(4)-3*u(5))/12);
s(N+2)=2*(exp(-1)-(25*u(N+1)-23*u(N)+13*u(N-1)-3*u(N-2))/12);
E = A\s;

x = x(2:N+1);
u = u(2:N+1);

ue = exp(-x.^2);
uea = uea(2:N+1);
plot(x,u,'o',x,ue,'+',x,uea,'*');
xlim([0 1])

figure
plot(x,ue-u,'x',x,uea-u,'s');
xlim([0 1])

%hold on
%plot(x,uea-ue)
%hold off

superr=max(abs(uea-u))
e = uea-u;



%figure 
%plot(x,u+e,x,ue);
%plot(x,u,'*')
%xlim([0 1])
%ubar = zeroes(N,1);
%for i = 1:10*N
%    if(
%   ubar(i)  
%end


E = E(2:N+1);




figure
plot(x,E,'v',x,e,'s','LineWidth',2)
xlim([0 1])
set(gca,'FontSize',18)
xlabel('x')
ylabel('error')
h=legend('Computed error average from transport equation', '$$\bar{u}$$-$$\tilde{u}$$, CV average of exact - computed');
set(h,'Interpreter','latex','fontsize',18)
errorerror=max(abs(E-e))
%figure
%max(abs(s(2:N+1)-f(2:N+1)))
%plot(x,s(2:N+1),x,f(2:N+1))

%max(abs(ue-q(2:N+1)))

%max(abs(uea-ue))




%%%
clear all
close all
N = 80;
s = 0.00001;%.001%1/N;
h0 = -s/(1-(1+s)^N);
h = zeros(N+2,1);
x = zeros(N+2,1);

for i = 1:N+2
    h(i) = h0*(1+s)^(i-2);
    if i == 1
        x(i) = -h(1)/2;
    elseif (i == 2)
        x(i) = h0/2;
    else
        x(i) = x(i-1)+(h(i-1)+h(i))/2;
    end
end
f = zeros(N+2,1);
for i = 1:N+2
    if i==1
       f(i) = 0; 
    else
    f(i) = (1/h(i))*(-2*(x(i)+h(i)/2)*exp(-(x(i)+h(i)/2)^2)+2*(x(i-1)+h(i-1)/2)*exp(-(x(i-1)+h(i-1)/2)^2));
    end
end

f(1)=2*1;
f(N+2)=2*exp(-1);


A = zeros(N+2,N+2);
for i = 1:N+2
   for j = 1:N+2
      if(i==j)
          if(i>1 && i < N+2)
         A(i,i) = (1/h(i))*(-1/((h(i+1)+h(i))/2)-1/((h(i-1)+h(i))/2)); 
          end
      end
      if(i+1==j)
          
         A(i,j) = (1/h(i))*(1/((h(i+1)+h(i))/2)); 
        
      end
      if(i-1==j)
          
         A(i,j)= (1/h(i))*(1/((h(i-1)+h(i))/2));
          
      end
   end
end


A(1,1) = 1;
A(1,2) = 1;
A(N+2,N+1)=1;
A(N+2,N+2)=1;

u=A\f;

uea = zeros(N+2,1);
for i = 1:N+2
    if(i > 1)
uea(i) = (1/h(i))*(sqrt(pi)*erf(x(i)+h(i)/2)/2-sqrt(pi)*erf(x(i-1)+h(i-1)/2)/2);
    end
end


E = zeros(N+2,1);
s = zeros(N+2,1);





for i = 1:N+2
    
    if ((i > 2) && (i < N+1))
        
     Mo = [ 1 1 1 1; 3*h(i+1)/2 1*h(i)/2 -1*h(i)/2 -3*h(i-1)/2; 7*h(i+1)^2/6 1*h(i)^2/6 1*h(i)^2/6 7*h(i-1)^2/6; 15*h(i+1)^3/24 1*h(i)^3/24 -1*h(i)^3/24 -15*h(i-1)^3/24];
     M = [ 1 1 1 1; 
         ((h(i+1)+h(i+2))^2-h(i+1)^2)/(2*h(i+2)) h(i+1)/2 -h(i)/2 ((-h(i))^2-(-h(i)-h(i-1))^2)/(2*h(i-1));
      ((h(i+1)+h(i+2))^3-h(i+1)^3)/(6*h(i+2)) h(i+1)^2/6  h(i)^2/6  ((-h(i))^3-(-h(i)-h(i-1))^3)/(6*h(i-1));
      ((h(i+1)+h(i+2))^4-h(i+1)^4)/(24*h(i+2)) h(i+1)^3/24 -h(i)^3/24 ((-h(i))^4-(-h(i)-h(i-1))^4)/(24*h(i-1))];
b = [ 0 1 0 0 ]';

%stencil weights
v = M\b;
q(1) = v(1);
q(5) = -v(4);
for j = 2:4
   q(j) = v(j)-v(j-1); 
end   
q = q/h(i);
        
        
    %s(i) =  (-u(i-2)+16*u(i-1)-30*u(i)+16*u(i+1)-u(i+2))/(12*h(i)^2);
    s(i) = (q(1)*u(i-2)+q(2)*u(i-1)+q(3)*u(i)+q(4)*u(i+1)+q(5)*u(i+2));
    end
    %if ((i > 1) && (i < N+2))
    %s(i) =  (u(i-1)-2*u(i)+u(i+1))/(h^2);
    %end
    
   % s(i) = (-2+4*((i-1.5)*h)^2)*(exp(-((i-1.5)*h)^2));
   
    %s(i) = s(i) - (1/h(i))*(-2*(x(i)+h(i)/2)*exp(-(x(i-1)+h(i-1)/2)^2)+2*(x(i-2)+h(i-2)/2)*exp(-(x(i-2)+h(i-2)/2)^2));
    if i > 1
    s(i) = s(i) -(1/h(i))*(-2*(x(i)+h(i)/2)*exp(-(x(i)+h(i)/2)^2)+2*(x(i-1)+h(i-1)/2)*exp(-(x(i-1)+h(i-1)/2)^2));
    end
   % end
end

j=2;

Mo = [ 1 1 1 1 1; 7*h(j+3)/2 5*h(j+2)/2 3*h(j+1)/2 1*h(j)/2 -1*h(j)/2 ;37*h(j+3)^2/6 19*h(j+2)^2/6 7*h(j+1)^2/6 1*h(j)^2/6 1*h(j)^2/6; 
    175*h(j+3)^3/24 65*h(j+2)^3/24 5*h(j+1)^3/8 1*h(j)^3/24 -1*h(j)^3/24; 781*h(j+3)^4/120 211*h(j+2)^4/120 31*h(j+1)^4/120 1*h(j)^4/120 1*h(j)^4/120 ];
M = [ 1 1 1 1 1;
    ((h(j+1)+h(j+2)+h(j+3)+h(j+4))^2-(h(j+1)+h(j+2)+h(j+3))^2)/(2*h(j+4))   ((h(j+1)+h(j+2)+h(j+3))^2-(h(j+1)+h(j+2))^2)/(2*h(j+3))   ((h(j+1)+h(j+2))^2-(h(j+1))^2)/(2*h(j+2))  h(j+1)/2 -h(j)/2;
    ((h(j+1)+h(j+2)+h(j+3)+h(j+4))^3-(h(j+1)+h(j+2)+h(j+3))^3)/(6*h(j+4))   ((h(j+1)+h(j+2)+h(j+3))^3-(h(j+1)+h(j+2))^3)/(6*h(j+3))   ((h(j+1)+h(j+2))^3-(h(j+1))^3)/(6*h(j+2))  h(j+1)^2/6 h(j)^2/6;
    ((h(j+1)+h(j+2)+h(j+3)+h(j+4))^4-(h(j+1)+h(j+2)+h(j+3))^4)/(24*h(j+4))  ((h(j+1)+h(j+2)+h(j+3))^4-(h(j+1)+h(j+2))^4)/(24*h(j+3))  ((h(j+1)+h(j+2))^4-(h(j+1))^4)/(24*h(j+2)) h(j+1)^3/24 -h(j)^3/24;
    ((h(j+1)+h(j+2)+h(j+3)+h(j+4))^5-(h(j+1)+h(j+2)+h(j+3))^5)/(120*h(j+4))  ((h(j+1)+h(j+2)+h(j+3))^5-(h(j+1)+h(j+2))^5)/(120*h(j+3))  ((h(j+1)+h(j+2))^5-(h(j+1))^5)/(120*h(j+2)) h(j+1)^4/120 h(j)^4/120];
    
b = [ 0 1 0 0 0 ]';
%stencil weights
v = M\b;
q(1) = v(1);
q(6) = -v(5);
for j = 2:5
   q(j) = v(j)-v(j-1); 
end
q = q(end:-1:1);
q = q/h(j);



%s(2)=s(2)+(10*u(1)-15*u(2)-4*u(3)+14*u(4)-6*u(5)+u(6))/(12*h(2)^2);
s(2)=s(2)+(q(1)*u(1)+q(2)*u(2)+q(3)*u(3)+q(4)*u(4)+q(5)*u(5)+q(6)*u(6));


j=N+1;
%Mo = [ 1 1 1 1 1; 7*h(j+3)/2 5*h(j+2)/2 3*h(j+1)/2 1*h(j)/2 -1*h(j)/2 ;37*h(j+3)^2/6 19*h(j+2)^2/6 7*h(j+1)^2/6 1*h(j)^2/6 1*h(j)^2/6; 
%    175*h(j+3)^3/24 65*h(j+2)^3/24 5*h(j+1)^3/8 1*h(j)^3/24 -1*h(j)^3/24; 781*h(j+3)^4/120 211*h(j+2)^4/120 31*h(j+1)^4/120 1*h(j)^4/120 1*h(j)^4/120 ];
M = [ 1 1 1 1 1;
    ((-h(j-1)-h(j-2)-h(j-3))^2-(-h(j-1)-h(j-2)-h(j-3)-h(j-4))^2)/(2*h(j-4))    ((-h(j-1)-h(j-2))^2-(-h(j-1)-h(j-2)-h(j-3))^2)/(2*h(j-3))   ((-h(j-1))^2-(-h(j-1)-h(j-2))^2)/(2*h(j-2))  -h(j-1)/2 h(j)/2;
    ((-h(j-1)-h(j-2)-h(j-3))^3-(-h(j-1)-h(j-2)-h(j-3)-h(j-4))^3)/(6*h(j-4))    ((-h(j-1)-h(j-2))^3-(-h(j-1)-h(j-2)-h(j-3))^3)/(6*h(j-3))   ((-h(j-1))^3-(-h(j-1)-h(j-2))^3)/(6*h(j-2))  h(j-1)^2/6 h(j)^2/6;
    ((-h(j-1)-h(j-2)-h(j-3))^4-(-h(j-1)-h(j-2)-h(j-3)-h(j-4))^4)/(24*h(j-4))   ((-h(j-1)-h(j-2))^4-(-h(j-1)-h(j-2)-h(j-3))^4)/(24*h(j-3))   ((-h(j-1))^4-(-h(j-1)-h(j-2))^4)/(24*h(j-2))  -h(j-1)^3/24 h(j)^3/24;
    ((-h(j-1)-h(j-2)-h(j-3))^5-(-h(j-1)-h(j-2)-h(j-3)-h(j-4))^5)/(120*h(j-4))   ((-h(j-1)-h(j-2))^5-(-h(j-1)-h(j-2)-h(j-3))^5)/(120*h(j-3))   ((-h(j-1))^5-(-h(j-1)-h(j-2))^5)/(120*h(j-2))  h(j-1)^4/120 h(j)^4/120];
b = [ 0 1 0 0 0 ]';
%stencil weights
v = M\b;
v = v(end:-1:1);
v = v/h(j);
q(1) = v(1);
q(6) = -v(5);
for j = 2:5
   q(j) = v(j)-v(j-1); 
end


%s(N+1)=s(N+1)+(10*u(N+2)-15*u(N+1)-4*u(N)+14*u(N-1)-6*u(N-2)+u(N-3))/(12*h(N+1)^2);
s(N+1)=s(N+1)+(q(1)*u(N+2)+q(2)*u(N+1)+q(3)*u(N)+q(4)*u(N-1)+q(5)*u(N-2)+q(6)*u(N-3));



s = -s;


j=0;
M = [ 1 1 1 1 ;
    ((h(j+1)+h(j+2)+h(j+3)+h(j+4))^2-(h(j+1)+h(j+2)+h(j+3))^2)/(2*h(j+4))   ((h(j+1)+h(j+2)+h(j+3))^2-(h(j+1)+h(j+2))^2)/(2*h(j+3))   ((h(j+1)+h(j+2))^2-(h(j+1))^2)/(2*h(j+2))  h(j+1)/2 ;
    ((h(j+1)+h(j+2)+h(j+3)+h(j+4))^3-(h(j+1)+h(j+2)+h(j+3))^3)/(6*h(j+4))   ((h(j+1)+h(j+2)+h(j+3))^3-(h(j+1)+h(j+2))^3)/(6*h(j+3))   ((h(j+1)+h(j+2))^3-(h(j+1))^3)/(6*h(j+2))  h(j+1)^2/6 ;
    ((h(j+1)+h(j+2)+h(j+3)+h(j+4))^4-(h(j+1)+h(j+2)+h(j+3))^4)/(24*h(j+4))  ((h(j+1)+h(j+2)+h(j+3))^4-(h(j+1)+h(j+2))^4)/(24*h(j+3))  ((h(j+1)+h(j+2))^4-(h(j+1))^4)/(24*h(j+2)) h(j+1)^3/24 ;];

b = [1 0 0 0]';
%stencil weights
v = M\b;
q = v;
q = q(end:-1:1);
q

%s(1) = 2*(1-(25*u(2)-23*u(3)+13*u(4)-3*u(5))/12);
s(1) = 2*(1-(q(1)*u(2)+q(2)*u(3)+q(3)*u(4)+q(4)*u(5)));

j=N+2;
M = [ 1 1 1 1 ;
    ((-h(j-1)-h(j-2)-h(j-3))^2-(-h(j-1)-h(j-2)-h(j-3)-h(j-4))^2)/(2*h(j-4))    ((-h(j-1)-h(j-2))^2-(-h(j-1)-h(j-2)-h(j-3))^2)/(2*h(j-3))   ((-h(j-1))^2-(-h(j-1)-h(j-2))^2)/(2*h(j-2))  -h(j-1)/2 ;
    ((-h(j-1)-h(j-2)-h(j-3))^3-(-h(j-1)-h(j-2)-h(j-3)-h(j-4))^3)/(6*h(j-4))    ((-h(j-1)-h(j-2))^3-(-h(j-1)-h(j-2)-h(j-3))^3)/(6*h(j-3))   ((-h(j-1))^3-(-h(j-1)-h(j-2))^3)/(6*h(j-2))  h(j-1)^2/6 ;
    ((-h(j-1)-h(j-2)-h(j-3))^4-(-h(j-1)-h(j-2)-h(j-3)-h(j-4))^4)/(24*h(j-4))   ((-h(j-1)-h(j-2))^4-(-h(j-1)-h(j-2)-h(j-3))^4)/(24*h(j-3))   ((-h(j-1))^4-(-h(j-1)-h(j-2))^4)/(24*h(j-2))  -h(j-1)^3/24 ];
b = [1 0 0 0]';
%stencil weights
v = M\b;
q = v;
q = q(end:-1:1);
q
%s(N+2)=2*(exp(-1)-(25*u(N+1)-23*u(N)+13*u(N-1)-3*u(N-2))/12);
s(N+2) = 2*(exp(-1)-(q(1)*u(N+1)+q(2)*u(N)+q(3)*u(N-1)+q(4)*u(N-2)));

E = A\s;

x = x(2:N+1);
h = h(2:N+1);
u=u(2:N+1);
plot(x,h,'*')
xlim([0 1])
uea = uea(2:N+1);
plot(x,uea,'*',x,u,'o')

superr=max(abs(uea-u))
e = uea-u;

E = E(2:N+1);


figure
plot(x,E,'v',x,e,'s','LineWidth',2)
xlim([0 1])
set(gca,'FontSize',18)
errorerror=max(abs(E-e))



% %%%% perturbed
% rng(1234);
% clear all 
% close all
% N=800;
% h0 = 1/N;
% X = zeros(N+1,1);
% for i = 1:N+1
%    X(i) = (i-1)*h0; 
%    if(i>1 && i < N+1)
%    X(i) = X(i) + (-1+rand*(2))*h0/3;
%    end
% end
% x = zeros(N+2,1);
% for i = 2:N+1
%     x(i) = (X(i-1)+X(i))/2;
% end
% x(1) = -x(2);
% x(N+2) = 1+(1-x(N+1));
% %plot(X,1,'*',x,1,'o')
% %xlim([0 1])
% h = zeros(N+2,1);
% for i = 2:N+1
%    h(i) = X(i)-X(i-1); 
% end
% 
% h(1) = h(2);
% h(N+2) = h(N+1);
% 
% f = zeros(N+2,1);
% for i = 1:N+2
%     if i==1
%        f(i) = 0; 
%     else
%     f(i) = (1/h(i))*(-2*(x(i)+h(i)/2)*exp(-(x(i)+h(i)/2)^2)+2*(x(i-1)+h(i-1)/2)*exp(-(x(i-1)+h(i-1)/2)^2));
%     end
% end
% 
% f(1)=2*1;
% f(N+2)=2*exp(-1);
% 
% 
% A = zeros(N+2,N+2);
% for i = 1:N+2
%    for j = 1:N+2
%       if(i==j)
%           if(i>1 && i < N+2)
%          A(i,i) = (1/h(i))*(-1/((h(i+1)+h(i))/2)-1/((h(i-1)+h(i))/2)); 
%           end
%       end
%       if(i+1==j)
%           
%          A(i,j) = (1/h(i))*(1/((h(i+1)+h(i))/2)); 
%         
%       end
%       if(i-1==j)
%           
%          A(i,j)= (1/h(i))*(1/((h(i-1)+h(i))/2));
%           
%       end
%    end
% end
% A(1,1) = 1;
% A(1,2) = 1;
% A(N+2,N+1)=1;
% A(N+2,N+2)=1;
% 
% u=A\f;
% 
% uea = zeros(N+2,1);
% for i = 1:N+2
%     if(i > 1)
% uea(i) = (1/h(i))*(sqrt(pi)*erf(x(i)+h(i)/2)/2-sqrt(pi)*erf(x(i-1)+h(i-1)/2)/2);
%     end
% end
% x = x(2:N+1);
% h = h(2:N+1);
% u=u(2:N+1);
% plot(x,h,'*')
% xlim([0 1])
% uea = uea(2:N+1);
% plot(x,uea,'*',x,u,'o')
% xlim([0 1])
% superr=max(abs(uea-u))
% e = uea-u;

