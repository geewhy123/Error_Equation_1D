
function [errerr2,x,cverr2,exacterr,ee ,te ] = errordriver( N,p,q,r ,unif,BCLeft,valLeft,BCRight,valRight,tlim,tord,physics,goal,jump,bchandle,NLflux)
ebcL = 0;
ebcR = 0;

if nargin == 13
    jump(1:3) = 0.2;
    bchandle = 'HC';
elseif nargin == 14
    bchandle = 'HC';
end

close all

if(p>0)
    rng(1234);
    
    %  g = randi(1000000);
    %  977219
    % rng(g)
    %   rng(972219);
    w = 0;
    h0 = 1/N;
    
    CFL = 0.4;
    k = CFL*h0;
    Q = tlim/(h0*CFL);
    if(abs(round(Q)-Q) > 1e-10)
        fprintf('not an integer number of steps, need to fix residual evaluation\n')
        error('1')
    end
    
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
    
    % load('tauN.mat')
    
    X = zeros(N+1,1);
    for i = 1:N+1%1:N+1
        X(i) = (i-1)*h0;
        if(i>1 && i < N+1)
            %    X(i) = X(i) + 0.1*randn*h0;
            X(i) = X(i) + unif*(-1+rand*(2))*h0/3;%0.001*sin(2*pi*X(i));%
            
            
        end
    end
    
    % X(2) = (2-1)*h0;
    % X(N) = (N-1)*h0;
    % X(N/2+1) = X(N/2+1)+(-1+rand*(2))*h0/3;
    
    x = zeros(N+2,1);
    for i = 2:N+1
        x(i) = (X(i-1)+X(i))/2;
        
        %     j = floor(i/2)+1;
        %     jj = mod(i,2);
        %     x(i) = x32(j)-h32(j)/4;
        %     if(jj == 1)
        %         x(i) = x(i) +h32(j)/2;
        %     end
    end
    
    x(1) = NaN;%0-(1-x(N+1));%-x(2);
    x(N+2) = NaN;%1+x(2);%1+(1-x(N+1));
    
    
    h = zeros(N+2,1);
    for i = 2:N+1
        h(i) = X(i)-X(i-1);
        
        %    j = floor(i/2)+1;
        %    h(i) = h32(j)/2;
    end
    
    h(1) = h(N+1);
    h(N+2) = h(2);
    
    % [x32 h32]
    % [x h]
    % plot(x,1,'*')
    % error('1')
    
    % global xx
    % xx = X;
    global dir
    dir = NaN*ones(N+2,1);
   
    if(strcmp(physics,'LinearSystem')==1)
        problem = pdelinearsystem(N,p,q,r,BCLeft,valLeft,BCRight,valRight,tlim,tord,physics,goal,x,h,k,0);
        problem.ebcL = ebcL;
        problem.ebcR = ebcR;
        problem.bchandle = bchandle;
    elseif(strcmp(physics,'EulerQ')==1)
        problem = pdeeuler(N,p,q,r,BCLeft,valLeft,BCRight,valRight,tlim,tord,physics,goal,x,h,k,0);
        problem.ebcL = ebcL;
        problem.ebcR = ebcR;
        problem.bchandle = bchandle;
        problem.NLfluxtype = NLflux;
    else
        problem = pde(N,p,q,r,BCLeft,valLeft,BCRight,valRight,tlim,tord,physics,goal,x,h,k,0);
        problem.jump = jump;
        problem.bchandle = bchandle;
    end
  
    problem.initializeexact();
    
    
    u0 = problem.initialSolution;
    ue = problem.exactSolution;
    length(x);
    length(ue);
    % figure
    %  plot(x,u0)
    
    %  plot(x,ue)
    %  assert(0)
    problem.computemoments();
    
    %u = ue;
    u=u0;
    %uu = zeros(N+2,1);
  
    %%%%
    % error('2')
    switch physics
        case {'Poisson', 'Advection'}
            problem.linearPhysics = true;
        case {'Burgers','BurgersVisc','EulerQ'}
            problem.linearPhysics = false;
        otherwise
            error('1')
            
    end
    
    if(problem.linearPhysics)
        problem.computeprimalpseudo();
        %       problem.jacobian = zeros(N+2,N+2);
        problem.primalJacobian = problem.computefluxjacobian(ue,'solution');
    end
   
    %%%%

    %   if((strcmp(physics,'Poisson')==1 && strcmp(goal,'SS')==1 && problem.bcLeftType == 'D' && problem.bcRightType == 'D' )||(strcmp(physics,'Advection')==1 && strcmp(goal,'SS')==1))
    if(strcmp(goal,'SS')==1 )
        fprintf('Implicit solve using Jacobian \n');
        %     [errerr2,x,cverr2,exacterr,ee,te  ]= problem.solvebyjacobian();
        if(strcmp(problem.physics,'EulerQ')~=1)
            [errerr2,x,cverr2,exacterr,ee,te  ]=problem.solvebyjacobianNL();
            return;
        end
        
        if(strcmp(problem.physics,'EulerQ')==1)
            [errerr2,x,cverr2,exacterr,ee,te  ]=solveeuler(problem);
            return;
        end
        return;
    end
 
    % global AD
    % AD = computepseudo(N,x,h,p,BCRight,BCLeft);
    problem.computeprimalpseudo();
    
    % problem.primalPI
  
    problem.curTime = 0;
    Z = problem.unstructuredrecon(ue,problem.pOrder,'solution');
 
    % figure
    if(strcmp(problem.physics,'LinearSystem')==1)
        er = problem.reconplot(Z(1:2,:),'solution');
        figure
        er = problem.reconplot(Z(3:4,:),'solution');
    end
    
    figure(1);clf;
    plt = plot(x,u0,'*');
    hold on
    grid on
    
    % Z
    % error('1')
    
    % hold on
    % plot(x,u0)
    f = problem.source;
    
    % figure
    % plot(x,ue)
    % error('1')
    
    Z;
    % % % % [tau]=reconfluxsoln(Z,f,h,N,p,physics,tlim,problem)
    [tau]=problem.computefluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,tlim,problem)
    
    
    if(strcmp(problem.physics,'LinearSystem')==1)
        problem.nUnk = 2;
    elseif(strcmp(problem.physics,'Advection')==1 ||strcmp(problem.physics,'Poisson')==1 || strcmp(problem.physics,'BurgersVisc')==1||strcmp(problem.physics,'Burgers')==1)
        
        problem.nUnk = 1;
    end
    
    
    d=1;
    u;
    % U =
    T = 1;
    klast = k;
    tlim
    problem.curTime = 0;
    
    
    Z = problem.unstructuredrecon(u,p,'solution');
    f = problem.computefluxintegral(Z,'solution');
    problem.Rall(:,1) = f;
    problem.Uall(:,1) = u;
    
    
    for j = 1:100000
        %     for m = 1:nUnk
        
        U(:,j,1:problem.nUnk) = u;
        tt = k*(j-1);
        
        
        if(tt+k > tlim)
            klast = tlim-tt;
            [tt k klast]
            if(klast < 1e-10)
                nSteps = j-1;
                T = (1:1:j-1)*k;
                T(end+1) = T(end)+klast;
                U(:,nSteps+1,1:problem.nUnk) = u;
                break
            end
            tt = tt +klast;
            
            problem.tStep = klast;
            %         problem.curTime = problem.curTime + klast;
            [uu,d] = problem.updatesolution(u);
            problem.curTime = problem.curTime + klast;
            u = uu;
            nSteps = j;
            T = (1:1:j-1)*k;
            T(end+1) = T(end)+klast;
            U(:,nSteps+1,1:problem.nUnk) = u;
            break
        end
        % if((max(d)*k<1e-15)  ||(tt>=tlim) )
        if((max(d)*k<1e-15))
            
            %         [uu,d] = problem.updatesolution(u);
            
            u = uu;
            max(d);
            tt;
            T = (1:1:j-1)*k;
            %  U(:,j+1) = u;
            nSteps = j-1;
            
            break
        end
        
        % d=0;      
        % problem.curTime = j*k;
        
        [uu,d] = problem.updatesolution(u);
        problem.curTime = j*k;
        
        
        % [uu,d] = update('solution',u,x,problem.source,k,h,N,p,tord,physics,NaN,NaN,problem);
        
        % if(j==20)
        %     uu
        % error('1')
        % end
        % [uu,d] = update('solution',u,x,f,k,h,N,p,tord,physics,NaN,NaN);  
        u = uu;
        
        
        set(plt,'ydata',u)
        drawnow
        
        %     end
    end
    % problem.Rall
    problem.Uall
    % U
    % error('1')
    % U
    fprintf('CFL = %e\n',CFL)
    fprintf('k = %e\n',k)
    fprintf('T_end = %e\n',tt)
    % nSteps
    % U(:,end-2:end)
    % pause
    u;
    [size(ue) size(u)]
    cverr1 = sum(abs(ue(2:N+1)-u(2:N+1)))/N;
    cverr2 = sqrt(sum((ue(2:N+1)-u(2:N+1)).^2)/N);
    cverrinf=max(abs(ue-u));
    fprintf('D.E.: [%e\t %e\t %e]\n\n',cverr1,cverr2,cverrinf);
    % size(U)
    
    
    u(1) = NaN;
    u(N+2) = NaN;
    plot(x,ue,'o')
    ylabel('u')
    figure
    plot(x,ue-u,'x')
    ylabel('ue-u')
    
end

tlim = T(end)

global dUdt
% if(klast > 1e-10 && klast < k)
dUdt = diffU(U,k,klast);
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


problem.convSoln = u;
if(q>0 && r > 0)
    
    clearvars -except u N p q r unif FI bta f cverr2 v k ue u0 tlim tord uo physics uder nSteps gsp U h x goal dUdt X problem Je tau

    
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
    
    problem.computerespseudo();%N,x,h,r);
    % problem.resPI
    % error('1')
    
    %global R
    R = zeros(N+2,nSteps+1);
    tt=0;
    
    
    % dUdt = dUdt*0;
    
    
    for j = 1:nSteps+1
        %    R(:,j) =computeres(U(:,j),x,h,N,f,r,physics,tt,gsp);
        R(:,j) = problem.computeres(U(:,j),tt,r);      
        tt = tt+k;        
    end
    
    
    if(exist('goal','var') && strcmp(goal,'FI')==1)
        assert(strcmp(physics,'Poisson')==1)
        for j = 1:nSteps+1
            R(:,j) = -FI;
        end
    end

    problem.residual = R;
    % problem.errorSource = R;
    % [-R(:,end) tau]
    % max(abs(-R(:,end)-tau))
    % plot(x,-R(:,end),x,tau)
    % R
    % error('1')
    
    U
    % R
    % error('1')
    
    plot(x,R(:,end))
    max(abs(R(:,end)));
    
    % dUdt
    %
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
    
    
    
    T=(0:1:nSteps)*k;
    for j = 2:N+1
        sp = spapi(6,T,R(j,:));
        Rsp(j) = sp;
        % spu = spapi(2,T,U(j,:));
        % spu = spapi(optknt(T,6),T,U(j,:));
        spu = pchip(T,U(j,:));
        Usp(j) = spu;
        % Rsp(j) = pchip(T,R(j,:));
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
    
    global M
    M = nSteps;
    fprintf('Error Equation\n')
    
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
   
    T = 1;
    s=1;
   
    figure(2);clf;
    plt2 = plot(x,u0,'*');
    hold on
    grid on
    axis([0 1 -1 1])
    % size(R)
    % R(:,end-2:end)

    problem.curTime = 0;
    for j = 1:nSteps+1
        
        
        E(:,j) = e;
        
        TT = k*(j-1);

        if(TT+k > tlim)
            %       TT=TT-k;
            [TT k tlim-TT]
            tlim-TT
            klast = tlim-TT;
            
            if(klast < 1e-10)
                
                nSteps = j-1;
                T = (1:1:j-1)*k;
                T(end+1) = T(end)+klast;
                
                E(:,nSteps+1,1:problem.nUnk) = e;
                break
            end
            
            problem.tStep = klast;
            problem.curTime = problem.curTime + klast;
            %             TT = TT +k;
            [ee,s] = problem.updateerror(e,TT,j);
            
            e = ee;
            nSteps = j;
            E(:,nSteps+1) = e;
            break
        end
  
        s=0;
        
        % % % [ee,s] = update('error',e,x,-R(:,j),k,h,N,q,tord,physics,TT,Rsp);
        % [e;TT]
        [ee,s] = problem.updateerror(e,TT,j);
        
     
        e = ee;        
        if(mod(j,100)==0)
            max(s)
        end
        
        set(plt2,'ydata',e)
        drawnow
     
        problem.curTime = j*k;
        
    end   
    exacterr = ue-u;

    exacterr = exacterr(2:N+1);
    x = x(2:N+1);
    figure
    plot(x,exacterr,'o-',x,ee(2:N+1),'*');
    
    
    ee = ee(2:N+1);
    ue=ue(2:N+1);
    
    
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


exacterr-ee;
te = NaN;

clear global
end
