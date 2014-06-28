   clear all
   close all
   N = 10;
    unif = 0;1/3;
   p = 2;
   q = 0;
   r = 0;
   BCLeft = 'D';
   BCRight = 'D';
   valLeft = 0;
   valRight = 0;
   tlim = 10;
   tord = 7;
   physics = 'Poisson';
   goal = 'SS';
    rng(1234);


w = 0;
h0 = 1/N;

CFL = 0.2;
k = CFL*h0;


    k = k*h0;

    k = 2*k;

     k = 1.5*k;




X = zeros(N+1,1);
for i = 1:N+1
   X(i) = (i-1)*h0; 
   if(i>1 && i < N+1)
%    X(i) = X(i) + 0.1*randn*h0;
   X(i) = X(i) + unif*(-1+rand*(2))*h0/3;%0.001*sin(2*pi*X(i));%
   end
end

x = zeros(N+2,1);
for i = 2:N+1
    x(i) = (X(i-1)+X(i))/2;
end

x(1) = NaN;%0-(1-x(N+1));%-x(2);
x(N+2) = NaN;%1+x(2);%1+(1-x(N+1));

h = zeros(N+2,1);
for i = 2:N+1
   h(i) = X(i)-X(i-1); 
end

h(1) = h(N+1);
h(N+2) = h(2);



problem = pde(N,p,q,r,BCLeft,valLeft,BCRight,valRight,tlim,tord,physics,goal,x,h,k,0);
problem.initializeexact();
f = problem.source;
 
u0 = problem.initialSolution;
ue = problem.exactSolution;
% plot(x,u0,x,ue)
% assert(0)
problem.computemoments();

%u = ue;
u=u0;
%uu = zeros(N+2,1);

problem.exactSolution



    fprintf('solving by Jacobian');
    [errerr2,x,cverr2,exacterr,ee,te  ]= problem.solvebyjacobian();

    
    u = problem.exactSolution - exacterr;
 
    g = zeros(N+2,1);
    ve = zeros(N+2,1);
    for i = 2:N+1
        xr = x(i)+h(i)/2;
        xl = x(i)-h(i)/2;
        g(i) = (1/h(i))*(pi^3/4)*(xr^2/2-xr^3/3-xl^2/2+xl^3/3);
        ve (i) = (1/h(i))*(-pi^3/48)*(xr^5/5-2*xr^4/4+xr^2/2-xl^5/5+2*xl^4/4-xl^2/2);
    end
    
    problem.exactSolution = ve;
    problem.source = g;
    
    [errerr2,x,cverr2,exacterr,ee,te  ]= problem.solvebyjacobian();
    
    v = ve-exacterr;
    
   order = 6; 
   problem.pOrder = order;
   problem.computemoments();
   problem.computeprimalpseudo();
    Zu = problem.unstructuredrecon(u,order,'solution');
    Zv = problem.unstructuredrecon(v,order,'solution');
    
%     error('1')
    J1 = 0;
    J2 = 0;
    
  
   c1 = 0.3478548451;
c2 = 0.6521451549;
c3 = 0.6521451549;
c4 = 0.3478548451;
x1= 0.8611363116;
x2 = 0.339981436;
x3 = -0.339981436;
x4= -0.8611363116;
G = @(x) (pi^3/4)*x*(1-x);
F = @(x) -pi^2*sin(pi*x);
    
    for i = 2:N+1
       xl = x(i)-h(i)/2;
        xr = x(i)+h(i)/2;      
 xx1 = ((xr-xl)/2)*x1+(xr+xl)/2;
 xx2 = ((xr-xl)/2)*x2+(xr+xl)/2;
 xx3 = ((xr-xl)/2)*x3+(xr+xl)/2;
 xx4 = ((xr-xl)/2)*x4+(xr+xl)/2;
%  error('2')
 U1 = 0;
 U2 = 0;
 U3 = 0;
 U4 = 0;
 V1 = 0;
 V2 = 0;
 V3 = 0;
 V4 = 0;
 for k = 1:order
 U1 = U1 + Zu(k,i)*(xx1-x(i))^(k-1); 
 U2 = U2 + Zu(k,i)*(xx2-x(i))^(k-1);
 U3 = U3 + Zu(k,i)*(xx3-x(i))^(k-1); 
 U4 = U4 + Zu(k,i)*(xx4-x(i))^(k-1);
 V1 = V1 + Zv(k,i)*(xx1-x(i))^(k-1); 
 V2 = V2 + Zv(k,i)*(xx2-x(i))^(k-1);
 V3 = V3 + Zv(k,i)*(xx3-x(i))^(k-1); 
 V4 = V4 + Zv(k,i)*(xx4-x(i))^(k-1);

 end
 [U1 U2 U3 U4];
 J1i = (c1*U1*G(xx1)+c2*U2*G(xx2)+c3*U3*G(xx3)+c4*U4*G(xx4))*(xr-xl)/2;
 J2i = (c1*V1*F(xx1)+c2*V2*F(xx2)+c3*V3*F(xx3)+c4*V4*F(xx4))*(xr-xl)/2;
 
 J1 = J1 + J1i;
        J2 = J2 + J2i;%v(i)*f(i)*h(i);
    
    end
    [J1 J2]