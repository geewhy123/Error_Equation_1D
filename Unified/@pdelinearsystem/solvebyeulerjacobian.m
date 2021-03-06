function  [errerr2,x,cverr2,exacterr,ee,te  ]  = solvebyeulerjacobian( obj)
%COMPUTEEULERJACOBIAN Summary of this function goes here
%   Detailed explanation goes here
    gam = obj.gamma;
    ue = obj.exactSolution;
    p = obj.pOrder;
    q = obj.qOrder;
    r = obj.rOrder;
    N = obj.nCells;
    x = obj.cellCentroids;

    V = obj.initialSolution;
    U = NaN*ones(3*N+2,1);
    obj.computeprimalpseudo();



    %%%%% higher order near bdy
    % [Z] = obj.unstructuredrecon(V,p,'solution');
    % 
    obj.hOrder = 0;%obj.pOrder;
    % m = obj.hOrder ;
    obj.computehigherpseudo();
    % 
    %       [Z3] = higherunstructuredreconeuler (obj,V(:,3),m,'solution',3);                          
    %                [Z1] = higherunstructuredreconeuler (obj,V(:,1),m,'solution',1); 
    %                 [Z2] = higherunstructuredreconeuler (obj,V(:,2),m,'solution',2);
    %                Zm = [Z1; Z2;Z3];
    % 
    % % error('1')



     %truncation error need exact sol
    teu = zeros(N+2,3);
    % tev = zeros(N+2,3);
    Ve = obj.exactSolutionV;
    Ue = obj.exactSolutionU;     
    if(obj.NLfluxtype==4)
           Ve = Ue;
    end

    % % % % higher
    [Z] = obj.unstructuredrecon(Ve,p,'solution');%ue,x,h,N,NaN,NaN,p);
    % % % % %  [Z3] = higherunstructuredreconeuler (obj,V(:,3),m,'solution',3);                          
    % % % % %                [Z1] = higherunstructuredreconeuler (obj,V(:,1),m,'solution',1); 
    % % % % %                 [Z2] = higherunstructuredreconeuler (obj,V(:,2),m,'solution',2);
    % % % % %                Z = [Z1; Z2;Z3];
    % % % % % % % % % % 


    %truncation error

    [tauu1 tauu2 tauu3]=obj.computeeulerfluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,tlim,obj)

    [tauu1 tauu2 tauu3];
    te1 = [ sum(abs(tauu1(2:N+1)))/N  sum(abs(tauu2(2:N+1)))/N sum(abs(tauu3(2:N+1)))/N ]

    teu = [tauu1 tauu2 tauu3];
 
    %truncation error need exact sol

  
    R = ones(3*N+2,1);
    total_time=0;
    u = ue;
    count = 0;
    c2 = 10;
    dtold = 0.01;
    Rold = R;
    
    if(obj.NLfluxtype==4)
       V = Ue; 
    end
    
    fprintf('Primal Equation\n')
    while(max(abs(R(2:3*N+1))) > 1e-13 )
        J = obj.computeeulerfluxjacobian(V,'solution');%,x,h,N,p);
    
        count = count +1;
         
        Rratio =norm(Rold(2:3*N+1),2)/norm(R(2:3*N+1),2); 
        dt = dtold*c2*Rratio;

%  spy(J)
%  J
%   error('1')


        K = J(2:3*N+1,2:3*N+1)+eye(3*N)/dt;

        [Z] = obj.unstructuredrecon(V,p,'solution');%u,x,h,N,NaN,NaN,p);
% % %%%%%
% % [Z] = obj.unstructuredrecon(V,p,'solution');%u,x,h,N,NaN,NaN,p);
% %   [Z3] = higherunstructuredreconeuler (obj,V(:,3),m,'solution',3);                          
% %   [Z1] = higherunstructuredreconeuler (obj,V(:,1),m,'solution',1); 
% %   [Z2] = higherunstructuredreconeuler (obj,V(:,2),m,'solution',2);
% %                Z = [Z1; Z2;Z3];
% % 
% % 
% % %%%%

        Rold = R;
        
        [phi1,phi2,phi3]=obj.computeeulerfluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)

% error('1')
        R(2:3:3*N-1) = phi1(2:N+1);
        R(3:3:3*N) = phi2(2:N+1);
        R(4:3:3*N+1) = phi3(2:N+1);
        u = NaN*ones(N+2,3);
        for j = 2:N+1
            [u(j,1),u(j,2),u(j,3)] = toconservedvars(V(j,1),V(j,2),V(j,3));
        end

        if(obj.NLfluxtype==4)
           u=V; 
        end

        del = K\-R(2:3*N+1);
    


        U(2:3:3*N-1) = u(2:N+1,1);
        U(3:3:3*N) = u(2:N+1,2);
        U(4:3:3*N+1) = u(2:N+1,3);

        UU = U(2:3*N+1) + del;
        U = NaN*ones(3*N+2,1);
        U(2:3*N+1) = UU;
          
        u(2:N+1,1) = U(2:3:3*N-1) ;
        u(2:N+1,2) = U(3:3:3*N) ;
        u(2:N+1,3) = U(4:3:3*N+1);

        for j = 2:N+1
            [V(j,1),V(j,2),V(j,3)] = toprimitivevars(u(j,1),u(j,2),u(j,3));
        end
        
        if(obj.NLfluxtype==4)
           V = u; 
        end
        
        total_time = total_time+dt;
          
        dtold = dt;
     
        
        fprintf('dt = %e , time = %e\n',  dt,total_time);
        fprintf('\t\t rho   Residual = %e \n' , max(abs(R(2:3:3*N-1))));
        fprintf('\t\t rho u Residual = %e \n' , max(abs(R(3:3:3*N-1))));
        fprintf('\t\t rho E Residual = %e \n' , max(abs(R(4:3:3*N-1))));
        
    end
  
    if(obj.NLfluxtype==4)
    for j = 2:N+1
             [V(j,1),V(j,2),V(j,3)] = toprimitivevars(u(j,1),u(j,2),u(j,3));
     end
    end
%       u
%     V
%     error('1')
    
    
    figure
    plot(x,V(:,1),'o',x,V(:,2),'v',x,V(:,3),'+')

    obj.convSoln = u;

    Ve = obj.exactSolutionV;    
   
    
    
    
    rho = obj.exactSolutionV(:,1);
    u = obj.exactSolutionV(:,2);
    P = obj.exactSolutionV(:,3);
    u1 = obj.exactSolutionU(:,1);
    u2 = obj.exactSolutionU(:,2);
    u3 = obj.exactSolutionU(:,3);
    uu = obj.convSoln;
    

    
    
    
    
    entropy = log(V(:,3)./(V(:,1).^gam));
    figure
    subplot(4,4,1)
    plot(x,V(:,1),'*',x,rho,'o')
    xlabel('$\rho$','Interpreter','Latex')
    subplot(4,4,2)
    plot(x,V(:,2),'*',x,u,'o')
    xlabel('u')
    subplot(4,4,3)
    plot(x,V(:,3),'*',x,P,'o')
    xlabel('P')
    subplot(4,4,4)
    plot(x,entropy,'*')
    xlabel('entropy')

    subplot(4,4,5)
    plot(x,rho-V(:,1),'x')
    subplot(4,4,6)
    plot(x,u-V(:,2),'x')
    subplot(4,4,7)
    plot(x,P-V(:,3),'x')
    subplot(4,4,8)
    plot(x,V(:,3)./V(:,1).^gam,'+')

    subplot(4,4,9)
    plot(x,uu(:,1),'*',x,u1,'o')
    xlabel('$\rho$','Interpreter','Latex')
    subplot(4,4,10)
    plot(x,uu(:,2),'*',x,u2,'o')
    xlabel('$\rho u$','Interpreter','Latex')
    subplot(4,4,11)
    plot(x,uu(:,3),'*',x,u3,'o')
    xlabel('$\rho E$','Interpreter','Latex')

    subplot(4,4,13)
    plot(x,u1-uu(:,1),'x')
    subplot(4,4,14)
    plot(x,u2-uu(:,2),'x')
    subplot(4,4,15)
    plot(x,u3-uu(:,3),'x')
    figure
    % plot(x,unew(:,1),x,unew(:,2),x,unew(:,3))
    plot(x,V(:,1),x,V(:,2),x,V(:,3))

    errerr2 = NaN;
    exacterrv = Ve-V;
    ee = NaN;
    cverr2 = [sqrt(sum((exacterrv(2:N+1,1)).^2)/N) sqrt(sum((exacterrv(2:N+1,2)).^2)/N) sqrt(sum((exacterrv(2:N+1,3)).^2)/N)];

    Ue = zeros(N+2,3);
    for i = 2:N+1
        [Ue(i,1),Ue(i,2),Ue(i,3)] = toconservedvars(Ve(i,1),Ve(i,2),Ve(i,3));
    end
    
    if(obj.NLfluxtype ==4)
        Ue = obj.exactSolutionU;
    end
    
    exacterr = Ue-obj.convSoln;
    cverru2 = [sqrt(sum((exacterr(2:N+1,1)).^2)/N) sqrt(sum((exacterr(2:N+1,2)).^2)/N) sqrt(sum((exacterr(2:N+1,3)).^2)/N)];
    

    fprintf ('\n\nt.e. rho   : %e   %e   %e\n' ,sum(abs(tauu1(2:N+1)))/N, sqrt(sum((tauu1(2:N+1)).^2)/N), max(abs(tauu1(2:N+1))));
    fprintf ('t.e. rho u : %e   %e   %e\n' ,sum(abs(tauu2(2:N+1)))/N, sqrt(sum((tauu2(2:N+1)).^2)/N), max(abs(tauu2(2:N+1))));
    fprintf ('t.e. rho E : %e   %e   %e\n\n' ,sum(abs(tauu3(2:N+1)))/N, sqrt(sum((tauu3(2:N+1)).^2)/N), max(abs(tauu3(2:N+1))));
    cverrv1 = exacterrv(2:N+1,1);
    cverrv2 = exacterrv(2:N+1,2);
    cverrv3 = exacterrv(2:N+1,3);
    
    fprintf ('d.e. rho : %e   %e   %e\n' ,sum(abs(cverrv1))/N, sqrt(sum((cverrv1).^2)/N), max(abs(cverrv1)));
    fprintf ('d.e.   u : %e   %e   %e\n' ,sum(abs(cverrv2))/N, sqrt(sum((cverrv2).^2)/N), max(abs(cverrv2)));
    fprintf ('d.e.   P : %e   %e   %e\n' ,sum(abs(cverrv3))/N, sqrt(sum((cverrv3).^2)/N), max(abs(cverrv3)));
        
    cverru1 = exacterr(2:N+1,1);
    cverru2 = exacterr(2:N+1,2);
    cverru3 = exacterr(2:N+1,3);
    
    fprintf ('\nd.e. rho : %e   %e   %e\n' ,sum(abs(cverru1))/N, sqrt(sum((cverru1).^2)/N), max(abs(cverru1)));
    fprintf ('d.e.   rho u : %e   %e   %e\n' ,sum(abs(cverru2))/N, sqrt(sum((cverru2).^2)/N), max(abs(cverru2)));
    fprintf ('d.e.   rho E : %e   %e   %e\n' ,sum(abs(cverru3))/N, sqrt(sum((cverru3).^2)/N), max(abs(cverru3)));
    te = teu;
    u=obj.convSoln;
    obj.convSolutionV = V;
    [Z] = obj.unstructuredrecon(V,p,'solution');
    obj.convVreconp = Z;

    obj.computeprimalleftright();

    [phi1,phi2,phi3]=obj.computeeulerfluxintegral(Z,'solution');


    %%%%%try to get BCs for error equation
    Pinf = 1;
    Tinf = 1;
    T0 = obj.T0;
    P0 = obj.P0;
    rhoinf = Pinf/Tinf;
    cinf = sqrt(gam*Pinf/rhoinf);
    uinf = cinf*sqrt((2/(gam-1))*(T0/Tinf-1));

    rhoi = 0;
    ui = 0;
    Pi = 0;

    h = obj.cellWidths;
    for k = 1:p
        rhoi = rhoi+ Z(k,2)*(-h(2)/2)^(k-1);
        ui   = ui+ Z(k+p,2)*(-h(2)/2)^(k-1);
        Pi   = Pi+ Z(k+2*p,2)*(-h(2)/2)^(k-1);
    end

    Vi= [rhoi ui Pi]';
    ci = sqrt(gam*Pi/rhoi);

    A = [cinf^2 0 -1; 0 -rhoinf*cinf -1; 0 rhoi*ci -1];
    b = [rhoinf-Pinf; -rhoinf*cinf*uinf-Pinf; rhoi*ci*ui-Pi];

    A\b;

    A\b-Vi;

    %%%%
    
    VV = obj.exactSolutionV;
    [Z] = obj.unstructuredrecon(VV,p,'solution');
 
    [Z(1,2)+Z(2,2)*-h(2)/2 Z(3,2)+Z(4,2)*-h(2)/2];
% error('1') 

 
    %%%%
    figure
    obj.reconplot(Z(1:p,:),'solution');
    figure
    obj.reconplot(Z(p+1:2*p,:),'solution');
    figure
    obj.reconplot(Z(2*p+1:3*p,:),'solution');
    save('tmp.mat','Z')
    % error('1')
% [obj.bcLeftVal(1) obj.bcLeftVal(2) obj.bcLeftVal(3)]

    Ze = obj.unstructuredrecon(Ve,p,'solution');
% Z
% 
Z = obj.unstructuredrecon(obj.convSolutionV,p,'solution');
(Ze(1,2)+Ze(2,2)*-0.1/2)-(Z(1,2)+Z(2,2)*-0.1/2)
(Ze(3,2)+Ze(4,2)*-0.1/2)-(Z(3,2)+Z(4,2)*-0.1/2)
% error('1')
% % Ze = Z;
% 0.957163318326954
% 0.961077533223941-(Ze(1,2)+Ze(2,2)*-0.1/2)
% 0.348597647123456
% 0.332088205777479-(Ze(3,2)+Ze(4,2)*-0.1/2)
% 0.97-(Ze(5,N+1)+Ze(6,N+1)*0.1/2)
% 

% 0.961077533223941-(Ze(1,2)+Ze(2,2)*-0.1/2+Ze(3,2)*(-0.1/2)^2+Ze(4,2)*(-0.1/2)^3)

% 0.332088205777479-(Ze(5,2)+Ze(6,2)*-0.1/2+Ze(7,2)*(-0.1/2)^2+Ze(8,2)*(-0.1/2)^3)

% error('1')
    
    
    %%%%residual

    if(q> 0 && r>0)

%         if(obj.NLfluxtype == 2)
%            obj.NLfluxtype = 4; 
%         end
        
        obj.computerespseudo();
    
        if(obj.NLfluxtype ==4)
           V = uu; 
        end
        Zr = obj.unstructuredrecon(V,r,'residual');

        figure
        obj.reconplot(Zr(1:r,:),'residual');
        figure
        obj.reconplot(Zr(r+1:2*r,:),'residual');
        figure
        obj.reconplot(Zr(2*r+1:3*r,:),'residual');

% [Zr] = obj.unstructuredrecon(u,r,'residual');
%   figure
%   obj.reconplot(Zr,'residual');
% error('1')
  
%     [left,right] = computeflux(Zr,h,N,r,physics,'residual',obj);
%     Rend= (right-left)./h-f;
        [R1, R2, R3] = obj.computeeulerfluxintegral(Zr,'residual');
 
        [R1 R2 R3];

        
        
                if(obj.NLfluxtype == 2 && strcmp(obj.bchandle,'HC')~= 1)
           obj.NLfluxtype = 4; 
        end


%         norm1R  = [sum(abs(R1(2:N+1)))/N sum(abs(R2(2:N+1)))/N sum(abs(R3(2:N+1)))/N]
%         error('1')

        if(obj.NLfluxtype == 4)
                   Ve = obj.exactSolutionU-uu; 
                   V= uu;
        end
        
        
        %%%error equation
        fprintf('\n\nError Equation\n')
        obj.computeerrorpseudo();
        %%%%
        [Zq] = obj.unstructuredrecon(Ve,q,'error');%ue,x,h,N,NaN,NaN,p);
        obj.reconexactsolutionV = Zq;
        %%%%

        Zu = obj.unstructuredrecon(V,obj.qOrder,'error');
        obj.convSolnRecon = Zu;
 
        if(obj.bcLeftType == 'D')
            obj.T0 = 0; 
            obj.P0 = 0;
        end
        if(obj.bcRightType == 'D')
            obj.Pb = 0;
        end

        exacterrv = obj.exactSolutionV-V;
        exacterru = obj.exactSolutionU-u;
        if(obj.NLfluxtype == 4)
            exacterrv = obj.exactSolutionV-obj.convSolutionV;
           
        end
        
        obj.exactSolutionV;
        UU =zeros(N+2,3);
        for z = 2:N+1
            [UU(z,1),UU(z,2),UU(z,3)] = toconservedvars(obj.exactSolutionV(z,1),obj.exactSolutionV(z,2),obj.exactSolutionV(z,3)); 
        end
        
        UU;
        obj.exactSolutionU;
% % error('1')

% obj.exactSolutionU
% error('1')

% figure
% plot(x,exacterr,'x')
%  error('1')

% % % need exact stuff
% % % %  obj.computeerrorpseudo();
% % % [Z] = obj.unstructuredrecon(exacterr,q,'error');%ue,x,h,N,NaN,NaN,p);
% % % 
% % % % Z
% % % % figure
% % % % obj.reconplot(Z,'error')
% % % % hold on
% % % % plot(x,exacterr)
% % % % error('1')


V = obj.exactSolutionV;
U = obj.exactSolutionU;




% load('tau4.mat')

        f = -[R1 R2 R3];

        obj.errorSource = f;teu;%f;%tau;
   
        figure
        subplot(3,1,1)
        plot(x,f(:,1),'+',x,teu(:,1),'^')
        subplot(3,1,2)
        plot(x,f(:,2),'+',x,teu(:,2),'^')
        subplot(3,1,3)
        plot(x,f(:,3),'+',x,teu(:,3),'^')
        legend('residual source','te source')
        [f teu];
        [ mean(abs(f(2:N+1,1))) mean(abs(f(2:N+1,2))) mean(abs(f(2:N+1,3)))];
%    error('1')
 
% % %  Je = obj.computefluxjacobian(exacterr,'error');
% % % % obj.errorRM
% % % % error('1')
% % % 
% % % 
% % % 
% % % % error('1')
% % % %  [tauE]=reconfluxsoln(Z,f,h,N,q,physics,tlim,obj)
% % % %  error('1')
% % %  [tauE]= obj.computefluxintegral(Z,'error')
% % % %   error('1')
       
        R = ones(3*N+2,1);
        total_time=0;
 
        e = exacterrv;
        
        ee = ones(3*N+2,1);
        ee(1)=NaN;
        ee(3*N+2) = NaN;
        count = 0;
        E = NaN*ones(3*N+2,1);
  
        dtold = 0.01;
        c2 = 10;

        Upe =zeros(N+2,3);
        Vpe = zeros(N+2,3);
        
        
        if(obj.NLfluxtype ~= 1)
           eu = exacterru;
           e = eu;
        end
        
%
        [Z] = obj.unstructuredrecon(e,q,'error');%u,x,h,N,NaN,NaN,p);
        
        [phi1 phi2 phi3]=obj.computeeulerfluxintegral(Z,'error');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)
%         phi2
%         error('1')
%      [Zvpe] = obj.unstructuredrecon(e+obj.convSolutionV,q,'error');%u,x,h,N,NaN,NaN,p);
%     [phi1 phi2 phi3]=obj.computeeulerfluxintegral(Zvpe,'error');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)

        obj.convVleft;
        obj.convVright;
    
        [phi1 phi2 phi3];
        figure
        plot(x,phi1,x,phi2,x,phi3)
%      error('1')
%
        R(2:3:3*N-1) = phi1(2:N+1);
        R(3:3:3*N) = phi2(2:N+1);
        R(4:3:3*N+1) = phi3(2:N+1);
        
        
        
        
%            delta = 1e-4*ones(size(e));
%         epd = e+delta;
%         [Z] = obj.unstructuredrecon(epd,q,'error');
%         [phi1 phi2 phi3]=obj.computeeulerfluxintegral(Z,'error');
%         R(2:3:3*N-1) = phi1(2:N+1);
%         R(3:3:3*N) = phi2(2:N+1);
%         R(4:3:3*N+1) = phi3(2:N+1);
%         
%         Je = obj.computeeulerfluxjacobian(e,'error');%,x,h,N,p);
%         [Z] = obj.unstructuredrecon(e,q,'error');
%         [phi1 phi2 phi3]=obj.computeeulerfluxintegral(Z,'error');
%         R0 = R;
%         R0(2:3:3*N-1) = phi1(2:N+1);
%         R0(3:3:3*N) = phi2(2:N+1);
%         R0(4:3:3*N+1) = phi3(2:N+1);
%         
%         Rdelta = R;
%    
%         Rdelta(2:3:3*N-1) = delta(2:N+1,1);
%         Rdelta(3:3:3*N) = delta(2:N+1,2);
%         Rdelta(4:3:3*N+1) = delta(2:N+1,3);
% 
%         R1 = R0+Je*Rdelta;
%         [R-R1]
%         error('1')

        
        while(max(abs(R(2:3*N+1))) > 1e-13  || total_time ==0)
     
% % % % if(obj.bcLeftType == 'D')
% % % %    obj.T0 = 0; 
% % % %    obj.P0 = 0;
% % % % end
% % % % if(obj.bcRightType == 'D')
% % % %     obj.Pb = 0;
% % % % end

            Je = obj.computeeulerfluxjacobian(e,'error');%,x,h,N,p);
% Je
% error('1')
% % % % % if(obj.bcLeftType == 'D')
% % % % %    obj.T0 = 1; 
% % % % %    obj.P0 = 1;
% % % % % end
% % % % % if(obj.bcRightType == 'D')
% % % % %     obj.Pb = 0.95;
% % % % % end

%      Jue = obj.computeeulerfluxjacobian(u+e,'error');%,x,h,N,p);
%      Ju = obj.computeeulerfluxjacobian(u,'error');%,x,h,N,p);
% Jue 
% Ju
%   Je
%        spy(Je)
%         error('1')
%      Je = Jue-Ju
% e
%      Je
%         error('1')
            count = count +1;

        
            Rratio =norm(Rold(2:3*N+1),2)/norm(R(2:3*N+1),2); 
    
            dt= dtold*c2*Rratio;

            K = Je(2:3*N+1,2:3*N+1)+eye(3*N)/dt;

            [Z] = obj.unstructuredrecon(e,q,'error');%u,x,h,N,NaN,NaN,p);

            Rold = R;
            [phi1 phi2 phi3]=obj.computeeulerfluxintegral(Z,'error');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)
   
            R(2:3:3*N-1) = phi1(2:N+1);
            R(3:3:3*N) = phi2(2:N+1);
            R(4:3:3*N+1) = phi3(2:N+1);
            
%     V
%     obj.convVreconp
%     [obj.convVleft obj.convVright]
% Z
%     [phi1 phi2 phi3]
%      figure
%      plot(x,phi1,x,phi2,x,phi3)
%     figure
%     plot(x,e)
%      error('1')

%%%%%new translation
if(obj.NLfluxtype == 1)
            Vpe = e+V;


            for j = 2:N+1
%     eu(j,1) = e(j,1);
%     eu(j,2) = e(j,2);
%     eu(j,3) = e(j,3);
% [eu(j,1),eu(j,2),eu(j,3)] = toconservedvars(e(j,1),e(j,2),e(j,3));
                [Upe(j,1),Upe(j,2),Upe(j,3)] = toconservedvars(Vpe(j,1),Vpe(j,2),Vpe(j,3));
            end

            eu = Upe-u;
% error('1')
            [Vpe Upe];
else
    eu = e;
end

%%%%%new translation
 
            del = K\-R(2:3*N+1);

            [Rratio max(abs(R(2:3*N+1)))];
    
            E(2:3:3*N-1) = eu(2:N+1,1);
            E(3:3:3*N) = eu(2:N+1,2);
            E(4:3:3*N+1) = eu(2:N+1,3);

            EE = E(2:3*N+1) + del;
            E = NaN*ones(N+2,1);
            E(2:3*N+1) = EE;
     
            eu(2:N+1,1) = E(2:3:3*N-1) ;
            eu(2:N+1,2) = E(3:3:3*N) ;
            eu(2:N+1,3) = E(4:3:3*N+1);

%%%%%%new trans
if(obj.NLfluxtype == 1)
            Upe = eu+u;

            for j = 2:N+1
%     e(j,1) = eu(j,1);
%     e(j,2) = eu(j,2);
%     e(j,3) = eu(j,3);
% [e(j,1),e(j,2),e(j,3)] = toprimitivevars(eu(j,1),eu(j,2),eu(j,3));
                [Vpe(j,1),Vpe(j,2),Vpe(j,3)] = toprimitivevars(Upe(j,1),Upe(j,2),Upe(j,3));
            end
         
            e = Vpe-V;
else 
    e = eu;
end
% error('1')
%%%%%%new trans


%      e
% exacterr

% error('1')
     
            total_time = total_time+dt;
            dtold = dt;
            
            
            fprintf('dt = %e , time = %e\n',  dt,total_time);
            fprintf('\t\t e_ rho   Residual = %e \n' , max(abs(R(2:3:3*N-1))));
            fprintf('\t\t e_ rho u Residual = %e \n' , max(abs(R(3:3:3*N-1))));
            fprintf('\t\t e_ rho E Residual = %e \n' , max(abs(R(4:3:3*N-1))));
            
        end


        if(obj.NLfluxtype ~=1)
        Upe = eu+obj.convSoln;
 
             for j = 2:N+1
     e(j,1) = eu(j,1);
     e(j,2) = eu(j,2);
     e(j,3) = eu(j,3);
 [e(j,1),e(j,2),e(j,3)] = toprimitivevars(eu(j,1),eu(j,2),eu(j,3));
                 [Vpe(j,1),Vpe(j,2),Vpe(j,3)] = toprimitivevars(Upe(j,1),Upe(j,2),Upe(j,3));
             end
          
             ee = Vpe-obj.convSolutionV;
        else
            ee = e;
        end
        

        w = exacterrv-ee;

% % % %  Je
% % % %  error('1')
% % % w = Je(2:N+1,2:N+1)\tauE(2:N+1)
% % % 
% % % % x1 = ones(N,1);
% % % %  null(Je(2:N+1,2:N+1),'r')
% % % % Je
% % % %   error('1')
% % % if(obj.bcLeftType=='P' && obj.bcRightType == 'P' && min(abs(eig(Je(2:N+1,2:N+1)))) < 1e-5)
% % % %   error('5')
% % % %   [Q,R] = qr(Je(2:N+1,2:N+1)') ;
% % % % w = Q*(R'\tauE(2:N+1)) ;
% % % % w = Je(2:N+1,2:N+1)'*Je(2:N+1,2:N+1)\Je(2:N+1,2:N+1)'*tauE(2:N+1)
% % % 
% % % % w = w-dot(w,x1')*x1*dot(w,w);
% % %     w = pinv(Je(2:N+1,2:N+1))*tauE(2:N+1)
% % %   end

        max(abs(w));
% error('1')
        figure
        plot(x,w)

% error('2')
 

% w(2:N+1) = w;
% w(1) = NaN;
% w(N+2) = NaN;
% cverr2 = sqrt(sum((ue(2:N+1)-u(2:N+1)).^2)/N);


% exacterr = exacterrv;
% exacterrv
% error('1')

% ee = exacterr - w;
        errerr2 = [sqrt(sum((exacterrv(2:N+1,1)-ee(2:N+1,1)).^2)/N) sqrt(sum((exacterrv(2:N+1,2)-ee(2:N+1,2)).^2)/N)  sqrt(sum((exacterrv(2:N+1,3)-ee(2:N+1,3)).^2)/N) ];
        w;
          figure
        subplot(3,2,1)
        plot(x,ee(:,1),'*',x,exacterrv(:,1),'o')
        xlabel('(computed)$\epsilon_\rho$','Interpreter','Latex');
        legend('estimate','"exact"')
        subplot(3,2,3)
        plot(x,ee(:,2),'*',x,exacterrv(:,2),'o')
        xlabel('(computed)$\epsilon_u$','Interpreter','Latex');
        subplot(3,2,5)
        plot(x,ee(:,3),'*',x,exacterrv(:,3),'o')
        xlabel('(computed)$\epsilon_P$','Interpreter','Latex');
        
        subplot(3,2,2)
        plot(x,eu(:,1),'*',x,exacterru(:,1),'o')
        xlabel('(translated)$\epsilon_\rho$','Interpreter','Latex');
        subplot(3,2,4)
        plot(x,eu(:,2),'*',x,exacterru(:,2),'o')
        xlabel('(translated)$\epsilon_{\rho u}$','Interpreter','Latex');
        subplot(3,2,6)
        plot(x,eu(:,3),'*',x,exacterru(:,3),'o')
        xlabel('(translated)$\epsilon_{\rho E}$','Interpreter','Latex');
        

        errerrv1 = exacterrv(2:N+1,1)-ee(2:N+1,1);
        errerrv2 = exacterrv(2:N+1,2)-ee(2:N+1,2);
        errerrv3 = exacterrv(2:N+1,3)-ee(2:N+1,3);
    
        fprintf ('d.e. e_ rho : %e   %e   %e\n' ,sum(abs(errerrv1))/N, sqrt(sum((errerrv1).^2)/N), max(abs(errerrv1)));
        fprintf ('d.e. e_  u : %e   %e   %e\n' ,sum(abs(errerrv2))/N, sqrt(sum((errerrv2).^2)/N), max(abs(errerrv2)));
        fprintf ('d.e. e_  P : %e   %e   %e\n' ,sum(abs(errerrv3))/N, sqrt(sum((errerrv3).^2)/N), max(abs(errerrv3)));
        
        errerru1 = exacterru(2:N+1,1)-eu(2:N+1,1);
        errerru2 = exacterru(2:N+1,2)-eu(2:N+1,2);
        errerru3 = exacterru(2:N+1,3)-eu(2:N+1,3);
        fprintf ('\nd.e. e_ rho : %e   %e   %e\n' ,sum(abs(errerru1))/N, sqrt(sum((errerru1).^2)/N), max(abs(errerru1)));
        fprintf ('d.e. e_  rho u : %e   %e   %e\n' ,sum(abs(errerru2))/N, sqrt(sum((errerru2).^2)/N), max(abs(errerru2)));
        fprintf ('d.e. e_  rho E : %e   %e   %e\n' ,sum(abs(errerru3))/N, sqrt(sum((errerru3).^2)/N), max(abs(errerru3)));
        
        
        
        errerru2 = [sqrt(sum((exacterru(2:N+1,1)-eu(2:N+1,1)).^2)/N) sqrt(sum((exacterru(2:N+1,2)-eu(2:N+1,2)).^2)/N)  sqrt(sum((exacterru(2:N+1,3)-eu(2:N+1,3)).^2)/N) ];
        errerr2 = max(abs(errerru2));
        save('euler244.mat','x','exacterru','eu')
    else
        errerr2 = NaN;
        exacterr = NaN;
        ee = NaN;
    
    end

end

