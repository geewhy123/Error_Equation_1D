
function [errerr2,x,cverr2,exacterr,ee ,te ] = errordriver( N,p,q,r ,unif,BCLeft,valLeft,BCRight,valRight,tlim,tord,physics,goal,jump,bchandle,NLflux,qRelin)
ebcL = 0;
ebcR = 0;

if nargin == 13
    jump(1:3) = 0.2;
    bchandle = 'HC';
    qRelin = 0;
elseif nargin == 14
    bchandle = 'HC';
    qRelin = 0;
end

close all

if(p>0)
    rng(1234);
    h0 = 1/N;
    
    CFL = 0.4;
%     CFL = 1;
% CFL = 0.2;

    k = CFL*h0;
% k = 1e-2;

    if(strcmp(physics,'Poisson')==1 || strcmp(physics,'BurgersVisc') == 1)
        
        if(strncmp(tord,'i',1)==1)
            k = k/4;
            
        else
            k = k*h0;
        end
        
        if( tlim > 2)
            k = 1.5*k;
        end
    end
    
%     k = 0.002;


    steps = tlim/k;
    if(abs(round(steps)-steps) > 1e-10)
        fprintf('not an integer number of steps, need to fix residual evaluation\n')
    end
    steps = round(steps);    
    X = zeros(N+1,1);
    for i = 1:N+1
        X(i) = (i-1)*h0;
        if(i>1 && i < N+1)
            X(i) = X(i) + unif*(-1+rand*(2))*h0/3;%0.001*sin(2*pi*X(i));
        end
    end
   
    x = zeros(N+2,1);
    for i = 2:N+1
        x(i) = (X(i-1)+X(i))/2;
    end
    
    x(1) = NaN;
    x(N+2) = NaN;
    h = zeros(N+2,1);
    for i = 2:N+1
        h(i) = X(i)-X(i-1);
    end

% w = 2;
% ff = @(s) 1+w*N*s-(1+s)^N;
% options = optimset(optimset('fsolve'), 'TolFun', 1.0e-16, 'TolX',1.0e-16,'MaxFunEvals',200);
% ss = fsolve(ff,0.5,options)
% assert(abs(ss) > 1e-12 )
% for j = 2:N+1
% h(j) = (1/(w*N))*(1+ss)^(j-2);
% end


% cc = 2;
% a = 1;
% w = ((1+cc/N^a)^N-1)/(cc*N^(1-a));
% for j = 2:N+1
% h(j) = (1/(w*N))*(1+cc/N^a)^(j-2);
% end

% w = 2;
% ff = @(s) 1-(2*w/N)*(1-1/(1+s)^(N/2))/(1-1/(1+s));
% options = optimset(optimset('fsolve'), 'TolFun', 1.0e-16, 'TolX',1.0e-16,'MaxFunEvals',200);
% ss = fsolve(ff,0.5,options)
% assert(abs(ss) > 1e-12 )
% for j = 2:N/2+1
% h(j) = (w/N)/(1+ss)^(N/2-j+1);
% end
% for j = N/2+2:N+1
% h(j) = (w/N)/(1+ss)^(j-N/2-2);
% end



for j = 2:N
X(j) = X(j-1)+h(j);
end

for j = 2:N+1
x(j) = 0.5*(X(j)+X(j-1));
end
% h
%  X
%  x
% cc/N
% error('1')


    h(1) = h(N+1);
    h(N+2) = h(2);    

% %nested subdivision
%     load('X10.mat')
%     ref = log(round(N/10))/log(2);
%     assert(mod(ref,1)==0);
%     
%     Xold = X10;
%     for jj = 1:ref
%         X = zeros(2*(length(Xold)-1)+1,1);
%         X(1:2:end) = Xold;
%         X(2:2:end-1) = (X(1:2:end-2)+X(3:2:end))/2;
%         Xold = X;
%     end

% % smooth 
% 
% s = 1/N;
% h0 = -s/(1-(1+s)^N)
% for i = 2:length(X)-1
%    X(i) = X(i-1)+h0*(1+s)^(i-2);
% %    X(i) = X(i) + +(-1+2*rand())*h0/3;
% end
% % X(N/2) = X(N/2)+(-1+2*rand())*h0/3;
%     for i = 2:N+1
%         h(i) = X(i)-X(i-1);
%     end
%     h(1) = h(N+1);
%     h(N+2) = h(2); 
%     for j = 2:N+1
%         x(j) = 0.5*(X(j)+X(j-1));
%     end
% %     X
% %     h
% %     x
% %     error('1')
    
    
    if(strcmp(physics,'LinearSystem')==1)
        problem = pdelinearsystem(N,p,q,r,BCLeft,valLeft,BCRight,valRight,tlim,tord,physics,goal,x,h,X,k,0);
        problem.ebcL = ebcL;
        problem.ebcR = ebcR;
        problem.bchandle = bchandle;
    elseif(strcmp(physics,'EulerQ')==1)
        problem = pdeeuler(N,p,q,r,BCLeft,valLeft,BCRight,valRight,tlim,tord,physics,goal,x,h,X,k,0);
        problem.ebcL = ebcL;
        problem.ebcR = ebcR;
        problem.bchandle = bchandle;
        problem.NLfluxtype = NLflux;
    else
        problem = pde(N,p,q,r,BCLeft,valLeft,BCRight,valRight,tlim,tord,physics,goal,x,h,X,k,0,qRelin);
        problem.jump = jump;
        problem.bchandle = bchandle;
         problem.NLError = NLflux;
    end
    
    if(strcmp(problem.physics,'BurgersVisc')==1 && strcmp(problem.goal,'TimeAccurate')==1)
%         problem.params.nu = 5e-3;
        problem.params.nu = 0.1;
    else
        problem.params.nu = 1;
    end
    
    problem.initializeexact();
    
    u0 = problem.initialSolution;
    ue = problem.exactSolution;
    
%     plot(x,u0,'o-',x,ue,'-+')
%     error('1')
    problem.computemoments();
    u=u0;
    switch physics
        case {'Poisson', 'Advection','Biharmonic'}
            problem.linearPhysics = true;
        case {'Burgers','BurgersVisc','EulerQ'}
            problem.linearPhysics = false;
        otherwise
            error('1')
            
    end
    
    if(problem.linearPhysics)
        problem.computeprimalpseudo();
        problem.primalJacobian = problem.computefluxjacobian(ue,'solution');
%         problem.primalJacobian
%         error('2')
    end
    
    if(strcmp(goal,'SS')==1 )
        fprintf('Implicit solve using Jacobian \n');
        if(strcmp(problem.physics,'EulerQ')~=1)
            [errerr2,x,cverr2,exacterr,ee,te  ]=problem.solvebyjacobianNL();
            if(q==0 && r ==0)
                errerr2 = cverr2;
                ee = NaN;
                exacterr = NaN;
            end
            return;
        end
        
        if(strcmp(problem.physics,'EulerQ')==1)
            [errerr2,x,cverr2,exacterr,ee,te  ]=solveeuler(problem);
            return;
        end
        return;
    end
    
    problem.computeprimalpseudo();
    problem.curTime = 0;
    Z = problem.unstructuredrecon(ue,problem.pOrder,'solution');
    if(strcmp(problem.physics,'LinearSystem')==1)
        er = problem.reconplot(Z(1:2,:),'solution');
        figure
        er = problem.reconplot(Z(3:4,:),'solution');
    end
    
    figure(1);clf;
    plt = plot(x,u0,'*');
    hold on
    grid on
    f = problem.source;
    
    [tau]=problem.computefluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,tlim,problem)
    
    
    if(strcmp(problem.physics,'LinearSystem')==1)
        problem.nUnk = 2;
    elseif(strcmp(problem.physics,'Advection')==1 ||strcmp(problem.physics,'Poisson')==1 || strcmp(problem.physics,'BurgersVisc')==1||strcmp(problem.physics,'Burgers')==1)
        
        problem.nUnk = 1;
    end
    
 
    
    klast = k;
    problem.curTime = 0;
    
    Z = problem.unstructuredrecon(u,p,'solution');
    f = problem.computefluxintegral(Z,'solution');
    problem.Rall(:,1) = f;
    problem.Uall(:,1) = u;
%     plot(x,problem.exactSolution,'o')
%     error('1')
% 
%   ZZ0 = problem.unstructuredrecon(problem.exactSolution(:,1),p,'solution');
%     ff0 = problem.computefluxintegral(ZZ0,'solution');
% ff0
% error('2')
    for j = 1:steps
        U(:,j,1:problem.nUnk) = u;
        
        [uu,d] = problem.updatesolution(u);
        
        problem.curTime = problem.curTime+k;%j*k;
        
        
        u = uu;
        
        set(plt,'ydata',u)
        drawnow
    end
    U(:,end+1,1:problem.nUnk) = u;

    % problem.Rall
    %     problem.Uall
    % U
    % error('1')
    % U
    fprintf('CFL = %e\n',CFL)
    fprintf('k = %e\n',k)
    fprintf('T_end = %e\n',problem.curTime)
    cverr1 = sum(abs(ue(2:N+1)-u(2:N+1)))/N;
    cverr2 = sqrt(sum((ue(2:N+1)-u(2:N+1)).^2)/N);
    cverrinf=max(abs(ue-u));
    fprintf('D.E.: [%e\t %e\t %e]\n\n',cverr1,cverr2,cverrinf);
    
    u(1) = NaN;
    u(N+2) = NaN;
    plot(x,ue,'o')
    ylabel('u')
    figure
    plot(x,ue-u,'x')
    ylabel('ue-u')
    
end

global dUdt
% if(klast > 1e-10 && klast < k)
dUdt = diffU(U,k,klast);

% dUdt = dUdt *0;
%  for j = 1:steps
%     dUdt(:,j) = (U(:,j+1)-U(:,j))/k;
%  end
% dUdt(:,steps+1) = (U(:,steps+1)-U(:,steps))/k;
 
% else
% dUdt = diffU(U,k);
% end
gsp = NaN;

% if not using FD, then use this to spline
% % % T=(0:1:nSteps)*k;
% % % for j = 2:N+1
% % % sp = spapi(6,T,U(j,:));
% % % gsp(j) = sp;
% % % %figure
% % % %fnplt(sp)
% % % %hold on
% % % %plot(T,U(3,:),'*')
% % % end

% save('test.mat','J','u0')
% error('1')
problem.convSoln = u;
[um un] = size(problem.Uall);
if(um == 0 || un == 0 )
   problem.Uall = u; 
end
if(q>0 && r > 0)
    
    clearvars -except u N p q r unif FI bta f cverr2 v k ue u0 tlim tord uo physics uder nSteps gsp U h x goal dUdt X problem Je tau steps
    

    figure
    hold on
    
    global KK
    KK = k;
    global xx
    xx = x;
    global TEND
    TEND = tlim;
    global UU
    UU = U;
   
    tt = 0;
    Rp = zeros(N+2,steps+1);%nSteps+1);
    for j = 1:steps+1%nSteps+1
        Z = problem.unstructuredrecon(U(:,j),p,'solution');
        Rp(:,j) =problem.computefluxintegral(Z,'solution');
        tt = tt+k;
    end
        
    problem.computerespseudo();%N,x,h,r);

    R = zeros(N+2,steps+1);%nSteps+1);
    tt=0;
    
%     dUdt
%     error('1')
%         dUdt = dUdt*0;
    % dUdt(2:N+1,:)=(problem.primalJacobian(2:N+1,2:N+1))*(U(2:N+1,:));

    for j = 1:steps+1%nSteps+1
       
        R(:,j) = problem.computeres(U(:,j),tt,r);
%  R(:,j) = problem.computeres(problem.exactSolution(:,1),tt,r);
        tt = tt+k;
    end
%     error('2')
%     U
%     R
%     error('1')
    
    
    if(exist('goal','var') && strcmp(goal,'FI')==1)
        assert(strcmp(physics,'Poisson')==1)
        for j = 1:nSteps+1
            error('1')
            R(:,j) = -FI;
        end
    end
    
    problem.residual = R;
    % problem.errorSource = R;
    % [-R(:,end) tau]
    % max(abs(-R(:,end)-tau))
    % plot(x,-R(:,end),x,tau)
    
    
    plot(x,R(:,end))
    ylabel('R')
    max(abs(R(:,end)));
    
    % t = (0:1:nSteps)*k;
    % Utexact = zeros(size(dUdt));
    % Utexact(1,:) = NaN;
    % Utexact(N+2,:) = NaN;
    % Uexact = Utexact;
    % for i = 2:N+1
    %     for j = 1:nSteps+1
    %           xl = x(i)-h(i)/2;
    %     xr = x(i)+h(i)/2;
    %         Uexact(i,j) = -(1/(2*pi))*(1/h(i))*(cos(2*pi*(xr+t(j)))-cos(2*pi*(xl+t(j))));
    %         Utexact(i,j) =           (1/h(i))*(sin(2*pi*(xr+t(j)))-sin(2*pi*(xl+t(j))));
    %     end
    % end
    %
    % Utexact
    % plot(x,dUdt(:,end)-Utexact(:,end))
    % max(abs(dUdt(:,end)-Utexact(:,end)))
    % % plot(x,Uexact(:,end))
    %
    %  error('1')
    
    
    
    T=(0:1:steps)*k;%nSteps)*k;
    for j = 2:N+1
        sp = spapi(2,T,R(j,:));
        Rsp(j) = sp;
        % spu = spapi(2,T,U(j,:));
        % spu = spapi(optknt(T,6),T,U(j,:));
        spu = pchip(T,U(j,:));
        Usp(j) = spu;
        % Rsp(j) = pchip(T,R(j,:));
%         fnval(k,Rsp(j))
%         fnval(1.5*k,Rsp(j))
%         fnval(2*k,Rsp(j))
%         error('1')
    end
    
    % figure
    % fnplt(Usp(2))
    % % U
    % % sum(U(2:N+1,:))
    % error('1')
    problem.Usp = Usp;
    
    e = zeros(N+2,1);
    ee =zeros(N+2,1);
    
    
    problem.Rsp =Rsp;
    
    
    figure
    Utall = problem.Uall;
    
    
    for i = 2:N+1
    Utall(i,:) = diffU(problem.Uall(i,:),k,k);
    end
    plot(x,Utall(:,end))
%     error('1')
    
    
    global M
    M = steps-1;%nSteps;
    fprintf('Error Equation\n')
    Emax = 0.1;%max(max(abs(problem.Uall-problem.exactSolutionAll)));
    % % %  AD = computepseudo(N,x,h,q);
    problem.computeerrorpseudo();
    
    Zu = problem.unstructuredrecon(problem.convSoln,problem.qOrder,'error');
    problem.convSolnRecon = Zu;
    
    
    
    if( problem.bcLeftType == 'D' && problem.bcRightType == 'D')
        problem.bcLeftVal = 0;
        problem.bcRightVal = 0;
    end
    
    if(problem.linearPhysics)
        timesbet = 0;
        for kk = 2:N+1
            val(kk,1:length(timesbet)) = fnval(timesbet,problem.Rsp(kk));
        end
        problem.errorSource = -1*val(:,1);
        problem.errorJacobian = problem.computefluxjacobian(e,'error');
    end
    
    
    
%     
%     RR = R;
%     J4=problem.errorJacobian;
%     save('test.mat','RR','J4')
    
    
    figure(2);clf;
    plt2 = plot(x,u0,'*');
    hold on
    grid on
    
    axis([0.1 1 -Emax Emax])
    
    problem.curTime = 0;
    for j = 1:steps%nSteps+1
        
        
        E(:,j) = e;
        s=0;
        [ee,s] = problem.updateerror(e,problem.curTime,j);%TT,j);
        e = ee;
        if(mod(j,100)==0)
            max(s)
        end
        
        set(plt2,'ydata',e)
        drawnow
        
        problem.curTime = problem.curTime+k;%j*k;
        
    end
    E(:,end+1,1:problem.nUnk) = e;
    
    
    exacterr = ue-u;
    
    exacterr = exacterr(2:N+1);
    x = x(2:N+1);
    figure
    plot(x,exacterr,'o-',x,ee(2:N+1),'*');
    
    
    ee = ee(2:N+1);
    ue=ue(2:N+1);
    
    problem.curTime
    errerr1 = sum(abs(exacterr-ee))/N;
    errerr2 = sqrt(sum((exacterr-ee).^2)/N);
    errerrinf=max(abs(exacterr-ee));
    fprintf('Error D.E.: [%e\t %e\t %e]\n',errerr1, errerr2, errerrinf);
    
    figure
    plot(x,exacterr-ee,'*-')
    ylabel('\epsilon - \epsilon_h')
    save('t','exacterr','ee','x')
    
else
    errerr2 = cverr2;
    ee = NaN;
    exacterr = NaN;
end
%     J = problem.primalJacobian;
%     J4 = problem.errorJacobian;
% %     save('jac.mat','J','U','E','R')
%     save('jac.mat','J')
%    [(ue(2:end-1)-u(2:end-1)) inv(J(2:end-1,2:end-1))*tau(2:end-1)]
%     
%     [norm(J(2:end-1,2:end-1),2) norm(inv(J(2:end-1,2:end-1)))]
%     error('1')
% error('1')
% exacterr
% ee
% E

exacterr-ee;
te = NaN;

% E=problem.exactSolutionAll-problem.Uall;
% save('tmp.mat','x','E')
% problem.source
% plot(x,problem.source(1,:))

% max(abs(ee))
clear global
end
