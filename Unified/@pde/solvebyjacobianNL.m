function [errerr2,x,cverr2,exacterr,ee,te  ] = solvebyjacobianNL( obj )
%SOLVEBYJACOBIAN Summary of this function goes here
%   Detailed explanation goes here

fprintf('\nPrimal Equation:\n')
ue = obj.exactSolution;
p = obj.pOrder;
q = obj.qOrder;
r = obj.rOrder;
h = obj.cellWidths;
N = obj.nCells;
x = obj.cellCentroids;


%%%%%%
% obj.hOrder = 0;
% obj.computehigherpseudo();
%%%%%%

obj.computeprimalpseudo();
[Z] = obj.unstructuredrecon(ue,p,'solution');%ue,x,h,N,NaN,NaN,p);

% obj.NLError = 'NLError';
obj.NLError = 'NewtonError';
disp(obj.NLError)

    
[tau]=obj.computefluxintegral(Z,'solution');
te1 = sum(abs(tau(2:N+1)))/N;
te2 = sqrt(sum((tau(2:N+1).^2)/N));
teinf = max(abs(tau(2:N+1)));

figure
plot(x,tau,'o')
te = tau;

f = obj.source;
% f = sourcefft(x,obj.source);
obj.source = f;

R = ones(N+2,1);
t=0;
u = ue;
count = 0;
c2 = 10;

total_time = 0;

Rold = R;
dtold = 1;
while(max(abs(R)) > 2e-11 )
    J = obj.computefluxjacobian(u,'solution');
    
%     [U,S,V] = svd(J(2:end-1,2:end-1))
%     error('5')
    
    count = count +1;
    if(count > 100)
        break;
    end
    Rratio =norm(Rold(2:N+1),2)/norm(R(2:N+1),2);
    dt = dtold*c2*Rratio;
    
    K = J(2:N+1,2:N+1)+eye(N)/dt;    
    [Z] = obj.unstructuredrecon(u,p,'solution');
    Rold = R;
    [R]=obj.computefluxintegral(Z,'solution');
    del = pinv(K)*-R(2:N+1);    
    uu = u(2:N+1) + del;
    u = NaN*ones(N+2,1);
    u(2:N+1) = uu;
    t = t+dt;
    
    dtold = dt;
    
    total_time = total_time+dt;
    fprintf('dt = %e , time = %e, Residual = %e \n', dt,total_time,max(abs(R)))
    
%     norm(inv(J(2:N+1,2:N+1))*N^2)
%     error('6')
    
end

vv = ue-u;
cverr1 = sum(abs(vv(2:N+1)))/N;
cverr2 =sqrt(sum((vv(2:N+1)).^2)/N);
cverrinf = max(abs(vv(2:N+1)));
fprintf('\nT.E.: [%e\t %e\t %e]\n',te1, te2, teinf);
fprintf('D.E.: [%e\t %e\t %e]\n',cverr1, cverr2, cverrinf);


plot(x,u,'*',x,ue,'o')

obj.convSoln = u;

% load('u244lin.mat')
% u = Upe;
% obj.convSoln = u;


if(strcmp(obj.goal,'SS')==1)
   obj.Uall = u; 
   obj.curTime = 0;
end
%h-truncation error
% Z = obj.unstructuredrecon(u,p,'solution')
% figure
% subplot(2,2,1)
% obj.reconplot(Z,'solution')
% ylabel('U_2')
% hold on
% obj.reconplot(Z,'average')
% ufine = zeros(2*N+2,1);
% for i = 2:N+1
%     for j = 0:p-1
%         ufine(2*(i-1)) = ufine(2*(i-1))+(2/h(i))*Z(j+1,i)*(-1)^j*(h(i)/2)^(j+1)/(j+1); 
%         ufine(2*(i-1)+1) = ufine(2*(i-1)+1)+(2/h(i))*Z(j+1,i)*(h(i)/2)^(j+1)/(j+1); 
%     end
% end
% u
% Z
%  obj.source
%  ufine
% % error('1')
% obj.meshrefinecoarsen('r');
% % obj.computerespseudo();
% Zfine = obj.unstructuredrecon(ufine,p,'solution')
% % error('1')
% % figure
% subplot(2,2,2)
% 
% obj.reconplot(Zfine,'average')
% axis([0 1 -3 3])
% hold on
% obj.reconplot(Zfine,'solution')
% % error('1')
% 
% % Rend = obj.computefluxintegral(Zfine,'residual');
% Rend = obj.computefluxintegral(Zfine,'solution');
%     subplot(2,2,4)
% plot(obj.cellCentroids,-Rend,'*')
% 
% Rend
% 
% Rhte = zeros(N+2,1);
% for i = 2:N+1
%    Rhte(i) = (Rend(2*(i-1))+Rend(2*(i-1)+1))/2;
% end
% Rhte

% obj.source
% % error('1')
% obj.meshrefinecoarsen('c')
% subplot(2,2,3)
% plot(x,-Rhte,'*')
% ylabel('h-\tau estimate')
%%%

if(q> 0 && r>0)
    
    fprintf('\n\n\nError Equation:\n')
    obj.computerespseudo();
    
    [Zr] = obj.unstructuredrecon(u,r,'residual');
    figure
    obj.reconplot(Zr,'residual');
    
    Rend = obj.computefluxintegral(Zr,'residual');
    
%     [Zre] = obj.unstructuredrecon(ue,r,'residual');
%     Rte = Rend;
%     Rte = Rte - obj.computefluxintegral(Zre,'residual');
% %     Rte = Rte +h(2)^5*rand(size(Rte));
% %     max(abs(Rte-Rend))
% [max(abs(Rte)) max(abs(Rend)) max(abs(Rte-Rend))]
%     error('1')
    obj.computeerrorpseudo();
    
    Zu = obj.unstructuredrecon(obj.convSoln,obj.qOrder,'error');
    obj.convSolnRecon = Zu;
    
    %%%%error equation
    uRL = 0;
    uRR = 0;
    for k = 1:p
        uRL = uRL + Z(k,2)*(-h(2)/2)^(k-1);
        uRR = uRR + Z(k,N+1)*(h(N+1)/2)^(k-1);
    end
    
    if(obj.bcLeftType == 'D')
        primalbcLeft = obj.bcLeftVal;
        obj.bcLeftVal = 0;uRL;(Z(1,2)+Z(2,2)*-h(2)/2); 0;
    end
    if(obj.bcRightType == 'D')
        primalbcRight = obj.bcRightVal;
        obj.bcRightVal = 0;uRR;(Z(1,N+1)+Z(2,N+1)*h(N+1)/2);0;
    end

    exacterr = ue-u;
    [Z] = obj.unstructuredrecon(exacterr,q,'error');%ue,x,h,N,NaN,NaN,p);
    
    f = -Rend;
 
    
%     f = -Rte;
%   f
%         Zr
%         error('1')
    
    
    
    mean(abs(f(2:N+1)-tau(2:N+1)));
    
%     figure
% %     subplot(1,2,1)
%     plot(x,-Rhte,'*',x,f,'o',x,tau,'^')
%     legend('R_{h-\tau}','R_{24}','\tau')
% %     subplot(1,2,2)
% %     plot(x,tau-f,'x')
    
tau0= tau;
% tau = tau+rand(N+2,1)*h(2)^1;
% load('te.mat')
[mean(abs(f(2:N+1)-tau(2:N+1))) max(abs(f(2:N+1)-tau(2:N+1)))]
plot(f)
hold on
plot(tau0)

% tau4 = tau;
% save('te.mat','tau4','-append')
% error('1')
    obj.errorSource = f;%tau4+f+rand(N+2,1)*h(2)^3;%f;%tau2-Rend;%tau2-Rend;f;%tau6-Rq;%f;%tau;
    
    [tauE]= obj.computefluxintegral(Z,'solution');
    tauE1 = sum(abs(tauE(2:N+1)))/N;
    tauE2 = sqrt(sum((tauE(2:N+1).^2)/N));
    tauEinf = max(abs(tauE(2:N+1)));
    
    R = ones(N+2,1);
    t=0;
    
    e = exacterr;
    ee = ones(N+2,1);
    ee(1)=NaN;
    ee(N+2) = NaN;
    count = 0;
    
    total_time = 0;
    
%     tau
% f
%     error('7')
e = 0*e;
    while(max(abs(R)) > 1e-11 )
        
        Je = obj.computefluxjacobian(e,'error');%,x,h,N,p);
        
        count = count +1;
        Rratio =norm(Rold(2:N+1),2)/norm(R(2:N+1),2);
        
        dt = dtold*c2*Rratio;
        K = Je(2:N+1,2:N+1)+eye(N)/dt;
        [Z] = obj.unstructuredrecon(e,q,'error');
        Rold = R;
        [R]=obj.computefluxintegral(Z,'error');
        
        
% % % % %         %
% % % % % %         norm(Je(2:N+1,2:N+1))
% % % % %          e(2:N+1) = Je(2:N+1,2:N+1)\-R(2:N+1);
% % % % %          e(N+2) = NaN;
% % % % % %          [e exacterr]
% % % % %          ee = e;
% % % % %          disp('newtonerror')
% % % % %          
% % % % % %          obj
% % % % %          %          error('1')
% % % % %          disp('re-linearize')
% % % % %          
% % % % % obj.bcLeftVal = primalbcLeft;
% % % % % obj.bcRightVal = primalbcRight;
% % % % %          obj.qOrder  = 6;
% % % % %          obj.rOrder = 6;
% % % % %           obj.convSoln = obj.convSoln+e;
% % % % % obj.computerespseudo();
% % % % % obj.computeerrorpseudo();
% % % % % obj
% % % % %           [Zs] = obj.unstructuredrecon(obj.convSoln,obj.qOrder,'error');
% % % % %          obj.convSolnRecon = Zs;
% % % % % 
% % % % % % error('1')
% % % % %          e = 0*e;
% % % % %          obj.bcLeftVal = 0;
% % % % %          obj.bcRightVal = 0;
% % % % %          Je6 = obj.computefluxjacobian(e,'error');
% % % % %          [norm(Je6(2:N+1,2:N+1)) norm(Je(2:N+1,2:N+1))]
% % % % %          [Z6] = obj.unstructuredrecon(e,obj.qOrder,'error');
% % % % %          [R6]=obj.computefluxintegral(Z6,'error');
% % % % %          e(2:N+1) = Je6(2:N+1,2:N+1)\-R6(2:N+1);
% % % % %          ee = e;
% % % % % %          full(Je)
% % % % %         break;
% % % % %         %
        
        
        del = pinv(K)*-R(2:N+1);
        
        max(abs(R(2:N+1)));
        
        ee(2:N+1) = e(2:N+1) + del;%*dt;
        e = NaN*ones(N+2,1);
        e = ee;
        t = t+dt;
        dtold = dt;
        
        total_time = total_time+dt;
        fprintf('dt = %e , time = %e, Residual = %e \n', dt,total_time,max(abs(R)))
    end
    
    
    w = exacterr-ee;

    ee = exacterr - w;
    
    errerr = exacterr-ee;
    errerr1 = sum(abs(errerr(2:N+1)))/N;
    errerr2 = sqrt(sum((errerr(2:N+1)).^2)/N) ;
    errerrinf = max(abs(errerr(2:N+1)));
    
    figure
    set(gca,'FontSize',25)
    h=plot(x,ee,'*',x,exacterr,'o','LineWidth',2)
    xlabel('x','FontSize',25)
%     ylabel('Discretization Error $$\epsilon$$','FontSize',30,'Interpreter','Latex')
%     legend('Estimate','Exact')
ylim([-0.04 0.03])
numberOfXTicks = 8;
xData = get(h,'XData');
set(gca,'Ytick',-0.04:0.01:0.03)
%set(gca,'YTickLabel',{'-0.04', '-0.03', '-0.02', '-0.01', '0', '0.01', '0.02', '0.03'})
    fprintf('\nError T.E.: [%e\t %e\t %e]\n',tauE1, tauE2, tauEinf);
    fprintf('Error D.E.: [%e\t %e\t %e]\n',errerr1, errerr2, errerrinf);
 
else
 
    errerr2 = NaN;
    exacterr = NaN;
    ee = NaN;
    
    
    
end
% 
    J = obj.primalJacobian;
% J4 = obj.errorJacobian;
% size(J)
% size(J4)
% %     save('jac.mat','J','U','E','R')
%     save('jac.mat','J4','-append')
% 
%    [norm((ue(2:end-1)-u(2:end-1))) norm(tau(2:end-1)) norm(inv(J(2:end-1,2:end-1))) norm(inv(Je(2:end-1,2:end-1))*J(2:end-1,2:end-1)) norm(inv(Je(2:end-1,2:end-1))*J(2:end-1,2:end-1))]
% J(2:end-1,2:end-1) = J(2:end-1,2:end-1)*Je(2:end-1,2:end-1);
%    [min(eig(J(2:end-1,2:end-1))) max(eig(J(2:end-1,2:end-1))) min(eig(J(2:end-1,2:end-1)*J(2:end-1,2:end-1)')) max(eig(J(2:end-1,2:end-1)*J(2:end-1,2:end-1)')) ] 
%    obj.primalJacobian
%     [norm(J(2:end-1,2:end-1),2)]
%     [norm(inv(J(2:end-1,2:end-1)))]
% load('jac.mat')
% [norm(J(2:end-1,2:end-1)*inv(J4(2:end-1,2:end-1)),2)]
% [norm(inv(J(2:end-1,2:end-1))*J4(2:end-1,2:end-1),2)]
% [norm(J(2:end-1,2:end-1)*inv(J4(2:end-1,2:end-1)),2)]
% [norm(inv(J4(2:end-1,2:end-1))*J(2:end-1,2:end-1),2)]
%     error('1')

% if q == 0
%  obj.convSoln
% else
% obj.convSoln+e
% end
if(obj.qRelin)
disp('re-linearize')
obj.bcLeftVal = primalbcLeft;
obj.bcRightVal = primalbcRight;
obj.qOrder  = 6;
obj.rOrder = 6;
obj.convSoln = obj.convSoln+e;

obj.computerespseudo();
[Z6] = obj.unstructuredrecon(obj.convSoln,obj.rOrder,'residual');
[R6]=obj.computefluxintegral(Z6,'residual');

obj.convSolnRecon = Z6;
obj.computeerrorpseudo();
[Zs] = obj.unstructuredrecon(obj.convSoln,obj.qOrder,'error');
% obj.convSolnRecon = Zs;

% error('1')
e = 0*e;
obj.bcLeftVal = 0;
obj.bcRightVal = 0;
Je6 = obj.computefluxjacobian(e,'error');
% Je6
% R6

% [norm(Je6(2:N+1,2:N+1)) norm(Je(2:N+1,2:N+1))]
e(2:N+1) = Je6(2:N+1,2:N+1)\-R6(2:N+1);



disp('')
disp('relin corr:')
% disp(max(abs(ue-(u+e))))
disp(max(abs(ue-(obj.convSoln+e))))



% figure
% plot(x,obj.convSoln+e)
% [ue obj.convSoln ]
%          full(Je)

%

 errerr = ue-(obj.convSoln+e);
    errerr2 = sqrt(sum((errerr(2:N+1)).^2)/N) ;
% Upe = obj.convSoln+e;
% save('u244lin.mat','Upe')

end

end

