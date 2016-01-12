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
    obj.computehigherpseudo();
    %%%%%%

    obj.computeprimalpseudo();


%     load('te.mat')
%     ue = u2040b;
    
    [Z] = obj.unstructuredrecon(ue,p,'solution');%ue,x,h,N,NaN,NaN,p);


 % % % % % % % Z
% % % % % % % % % error('2')
% % % % % % % %   [er]=obj.reconplot(Z,'solution')%reconplot(x,h,N,p,Z)
% % % % % % % % %   error('1')
% % % % % % % %   
% % % % % % % %   
% % % % % % % % f = obj.source;
% % % % % % % %  [tau]=obj.computefluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,tlim,obj)
% % % % % % % %  max(abs(tau))
% % % % % % % %  
% % % % % % % % %  sqrt(sum((tau(2:N+1)).^2)/N)
% % % % % % % % te = tau;
% % % % % % % % % 
% % % % % % % % % tau2 = tau;
% % % % % % % % % save('tau.mat','tau2','-append')
% % % % % % % % % error('1')
% % % % % % % % 
% % % % % % % % 
% % % % % % % % [sum(abs(tau(2:N+1)))/N sqrt(sum((tau(2:N+1).^2)/N)) max(abs(tau(2:N+1)))]
% % % % % % % % 
% % % % % % % % % errerr2= NaN;
% % % % % % % % % cverr2 = NaN;
% % % % % % % % % exacterr = NaN;
% % % % % % % % % ee = NaN;
% % % % % % % % % return;
% % % % % % % % figure
% % % % % % % % plot(x,tau,'o')
% % % % % % % % te=  sum(abs(tau(2:N+1)))/N;


obj.NLError = 'NLError';



    obj.hOrder = obj.pOrder;

    refinecells = [];%[2 3 4 N-1 N N+1];%2 3 4 5 N-2 N-1 N N+1];
    obj.refinecells = refinecells;
    if(obj.hOrder > 0)
        obj.computehigherpseudo();
        [Zh] = obj.unstructuredrecon(ue,obj.hOrder,'solution'); 
        Znew = zeros(max(obj.hOrder,p),N+2);
        Znew(1:p,:) = Z;
        
        if(obj.hOrder >= p)
            for ii = 1:length(refinecells)
                Znew(:,refinecells(ii)) = Zh(:,refinecells(ii));
            end
        else
            for ii = 1:length(refinecells)
                Znew(1:obj.hOrder,refinecells(ii)) = Zh(1:obj.hOrder,refinecells(ii));
                Znew(obj.hOrder+1:end,refinecells(ii)) = zeros(p-obj.hOrder,1);
            end
        end
        
        figure
        obj.reconplot(Z(1:p,:),'solution');
%         error('1')

%      c = -[1 -2 1]/h(3)^2;
        c = -[1 0 -1]/(2*h(3));
        d = [Znew(1,2) Znew(1,3) Znew(1,4)];
        cc = -[3 -4 1]/(2*h(2));
%     Znew(2,2) = dot(c,d);

%     Znew(2,2) = dot(cc,d);
        ccc = [-43/24 69/24 -33/24 7/24]/(h(2));
        ddd = [Znew(1,2) Znew(1,3) Znew(1,4) Znew(1,5)];
% % % Znew(2,2) = dot(ccc,ddd);

% Znew(:,3)
% error('1')
   
% Znew(2,3) = dot(c,d);%2.987132;dot(c,d);

        Zgradexact = Znew;
        Zstruct = Znew;
        for kk = 2:N+1
            Zgradexact(2,kk) =pi*cos(pi*x(kk));
            if(kk == 2)
                Zstruct(2,kk) = dot(cc,Znew(1,kk:kk+2));
            elseif(kk==N+1)
                Zstruct(2,kk) = dot(-cc(end:-1:1),Znew(1,kk-2:kk));
            else
                Zstruct(2,kk) = dot(c,Znew(1,kk-1:kk+1));
            end
        end
% Znew(2,:) = Zgradexact(2,:);
% Znew(2,2) = dot(ccc,Znew(1,2:5));
% Znew(2,2) = dot(cc,Znew(1,2:4));

        [Znew(2,4)  dot(c,Znew(1,3:5)) dot([5/48 -34/48 0 34/48 -5/48]/h(4),Znew(1,2:6)) pi*cos(pi*x(4))]    ;
%   [Znew(2,4)  dot(ccc,Znew(1,4:7)) pi*cos(pi*x(4))]    
%     Znew(2,2) = pi*cos(pi*h(2));
        err = abs(Znew(2,2)-pi*cos(pi*h(2)));
       
        Znew;
        Zgradexact;
        Zstruct;
%     Znew = Zstruct
    
%     error('1')
    end
    
    [tau]=obj.computefluxintegral(Znew,'solution');
    te1 = sum(abs(tau(2:N+1)))/N;
    te2 = sqrt(sum((tau(2:N+1).^2)/N));
    teinf = max(abs(tau(2:N+1)));

    figure
    plot(x,tau,'o')
%     tau;
% %     tau2040b = tau;
% %     save('te.mat','tau2040b','-append')
% dot(tau(2:N+1),h(2:N+1))
% dot(tau(N/2:N/2+3),h(N/2:N/2+3))
% tau(N/2+2)+tau(N/2+1)

%  x
%     error('1')

% %     error('1')
%     errerr2 =max(abs(tau(2:N+1)))
%     cverr2 = 0;
%     exacterr = 0;
%     ee = 0;
%     te = 0;
%     obj.jump
%     return;


%     hold on
%     a = tau;
%     b = a.*abs(a)/max(abs(a));
% plot(x,b)
%     
% mean(abs(tau-b))    
%     error('1')
    figure
    % obj.pOrder = obj.hOrder;
    obj.reconplot(Znew,'solution');
    % obj.reconplot(Z,'solution')
    % error('1') 
    Znew;
    te = tau;

    %   [er]=obj.reconplot(Z,'solution')%x,h,N,p,Z);
    %   error('1')
    f = obj.source;
    %  [tau]=obj.computefluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,tlim,obj)
%  save('tau.mat','x','tau')
%  error('1')
f = sourcefft(x,obj.source);
obj.source = f;

figure

    %   tau
% %gradient accuracy
%      Z
%         for j = 1:N
%             FR(j) = 0;
%             FL(j) = 0;
%             for k = 1:p-1
%            FR(j) = FR(j) + k*Z(k+1,j+1)*(h(j+1)/2)^(k-1);%Z(1,j+1)+Z(2,j+1)*(h(j+1)/2); 
%            FL(j) = FL(j) + k*Z(k+1,j+1)*(-h(j+1)/2)^(k-1);%Z(1,j+1)+Z(2,j+1)*(-h(j+1)/2); 
%             end
%         end
%         for j = 2:N
%            F(j) = (FL(j)+FR(j-1))/2; 
%            Fe(j) = pi*cos(pi*(x(j)+h(j)/2));
%         end
%         F(1) = FL(1);
%         F(N+1) = FR(N);
%         Fe(1) = pi*cos(pi*0);
%         Fe(N+1) = pi*cos(pi*1);
%         FL
%         FR
%         F
%         Fe
%         mean(abs(Fe-F))

    te = tau;

%       tau2 = tau;
%       save('tauB.mat','tau2')%,'-append')
%       error('1')

    % errerr2= NaN;
    % cverr2 = NaN;
    % exacterr = NaN;
    % ee = NaN;
    % return;

    %   error('1')
    
    R = ones(N+2,1);
    t=0;
    u = ue;
    count = 0;
    c2 = 10;
 
    total_time = 0;
 
    Rold = R;
    dtold = 1;
    while(max(abs(R)) > 2e-11 )
        J = obj.computefluxjacobian(u,'solution');%,x,h,N,p);
%         norm(J)
%         error('1')

%      q = eig(J);
%      req = real(q);
%      max(abs(req))
%      q
%      figure
%      plot(real(q),imag(q),'o')
% eig(J)
% error('1')
%      error('1')
        count = count +1;
        if(count > 100)
            break;
        end
        Rratio =norm(Rold(2:N+1),2)/norm(R(2:N+1),2); 
        dt = dtold*c2*Rratio;

        K = J(2:N+1,2:N+1)+eye(N)/dt;
% K
% error('1')
%  K = (K+K')/2;

        [Z] = obj.unstructuredrecon(u,p,'solution');%u,x,h,N,NaN,NaN,p);

       
%         error('1')
%  [er]=reconplot(x,h,N,p,Z);
        Rold = R;
        [R]=obj.computefluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)
%         del = (K'*K)\(K'*-R(2:N+1));%pinv(K)*-R(2:N+1);%K\-R(2:N+1);


%         del = K\-R(2:N+1);
del = pinv(K)*-R(2:N+1);
        
%     if(mod(count,100)==0)
    
%     end

        uu = u(2:N+1) + del;%*dt;
        u = NaN*ones(N+2,1);
        u(2:N+1) = uu;
        t = t+dt;
          
        dtold = dt;
     
        total_time = total_time+dt;
        fprintf('dt = %e , time = %e, Residual = %e \n', dt,total_time,max(abs(R)))
     
    end

    vv = ue-u;
    cverr1 = sum(abs(vv(2:N+1)))/N;
    cverr2 =sqrt(sum((vv(2:N+1)).^2)/N);
    cverrinf = max(abs(vv(2:N+1)));
    fprintf('\nT.E.: [%e\t %e\t %e]\n',te1, te2, teinf);
    fprintf('D.E.: [%e\t %e\t %e]\n',cverr1, cverr2, cverrinf);
  
    subplot(211)
    plot(x,u,'*',x,ue,'o')
    
%     error('1')
%     cverr2 = te2;
%     errerr2 = te2;
%     exacterr = te2;
%     ee = te2;
%     return ;
%    error('1')
%   u-ue

%  v=J(2:N+1,2:N+1)\tau(2:N+1)
% 
%   if(obj.bcLeftType=='P' && obj.bcRightType == 'P' && min(abs(eig(J(2:N+1,2:N+1)))) < 1e-5)
%     v = pinv(J(2:N+1,2:N+1))*tau(2:N+1);
%   end
% %  max(abs(v-ue(2:N+1)))
% max(abs(v))
% %  error('1')
%  figure
%  plot(x,u-ue,x(2:N+1),v)
%  
%  v(2:N+1) = v;
%  v(1) = NaN;
%  v(N+2) = NaN;
%  
%  u = ue-v
%  plot(x,u)
%  v
% %  ue-v
% % cverr1 = sum(abs(v(2:N+1)))/N
% %  cverr2 = sqrt(sum((v(2:N+1)).^2)/N)

 obj.convSoln = u;
%  u
% vv
%  obj.convSoln
%  error('1')
% save('up.mat','u')
% error('1')


% % % % accuracy
% % % % pp=obj.pOrder;
% % % % obj.pOrder = obj.qOrder;
% % % % p = obj.pOrder;
% % % % s = f;
% % % % obj.hOrder = 0;obj.pOrder;
% % % % obj.computeprimalpseudo();
% % % %  [Z] = obj.unstructuredrecon(ue,p,'solution');%ue,x,h,N,NaN,NaN,p);
% % % %  [tauq]=obj.computefluxintegral(Z,'solution');
% % % % 
% % % %  tauq
% % % %  obj.source = obj.source*0;
% % % %   [Z] = obj.unstructuredrecon(u,q,'solution');%ue,x,h,N,NaN,NaN,p);
% % % %  FIq = obj.computefluxintegral(Z,'solution');
% % % %  FIq
% % % % %  obj.computerespseudo();
% % % % obj.pOrder = r;
% % % %  obj.computeprimalpseudo();
% % % %   [Z] = obj.unstructuredrecon(u,r,'solution');%ue,x,h,N,NaN,NaN,p);
% % % %  FIr = obj.computefluxintegral(Z,'solution'); 
% % % %   FIr
% % % %  
% % % % 
% % % % % mean(abs(tauq(2:N+1)))
% % % % g=tauq-(FIq-FIr)
% % % % 
% % % % J = obj.computefluxjacobian(ue,'solution');%,x,h,N,p);
% % % % J(2:N+1,2:N+1)\(tauq(2:N+1)-FIq(2:N+1)+s(2:N+1))
% % % % J(2:N+1,2:N+1)\(-FIr(2:N+1)+s(2:N+1))
% % % % ue-u
% % % % % -FIr+s
% % % % mean(abs(g(2:N+1)))
% % % % 
% % % % b = tauq+s-FIq
% % % % mean(abs(b(2:N+1)))
% % % % c = FIr-s
% % % % mean(abs(c(2:N+1)))
% % % % 
% % % %  error('1')
% % % % obj.pOrder = pp;
% % % % p = pp;
% % % % obj.source = s;
% % % % 
%   J = obj.computefluxjacobian(0.5*(ue+u),'solution');%,x,h,N,p);
% J(2:N+1,2:N+1)*vv(2:N+1)
% tau

%     obj.postprocesscn(0.5);

%%%%residual
% J(2:N+1,2:N+1)%*u(2:N+1) - obj.source(2:N+1)
% obj.computeerrorpseudo();
% obj.errorSource = zeros(N+2,1);
% K = obj.computefluxjacobian(u,'error');
% J-K
% error('1')
% ue-u
% figure
%     plot(x,tau,'*')
%      save('tmp.mat','x','tau')

% v=J(2:N+1,2:N+1)*(ue(2:N+1)-u(2:N+1))
% 
% K = obj.computefluxjacobian(ue,'solution');
% w=K(2:N+1,2:N+1)*(ue(2:N+1)-u(2:N+1))
% (v+w)/2
% tau
% error('1')

%  u = u-mean(u(2:end-1));
%  hold on 
%  plot(x,u,'v')
% % % 



% hold on
% plot(x,obj.exactNoisySolution,'^')
% subplot(212)
% plot(x,ue-u,'x')
% ylabel('u_singleexact-u_smoothed')
% saveas(gca,'smootherror64b.eps','epsc')
% % error('1')
% figure
% plot(x,obj.source,'*',x,-4*pi^2*u,'+')
% mean(obj.source(2:end-1))
% 
% 
% mean(u(2:end-1))
% 
% 
% % error('1')
% x = x(2:end-1);
% u = u(2:end-1);
% L = length(x);
% NFFT = 2^nextpow2(L);
% Fs = 1*(length(x));
% f = Fs/2*linspace(0,1,NFFT/2+1);
% Y = fft(u,NFFT);
% figure
% subplot(223)
% stem(f,(2/L)*abs(Y(1:NFFT/2+1)));
% xlabel('fft of conv. soln')
% subplot(224)
% stem(f,abs(Y(1:NFFT/2+1)).*(2/L).*(2*pi*f').^2) 
% xlabel('fft of conv. soln * (2\pi f)^2')
% length(f')
% length(Y(1:NFFT/2+1))
% ylim([0 max(abs(Y(1:NFFT/2+1)).*(2/L).*(2*pi*f').^2)])
% (2)*abs(Y(1:NFFT/2+1)).*(2*pi*f').^2;
% 
% 
% tau64=tau;
% x64 = x;
% h64 = h;
% save('tauN.mat','x64','tau64','h64','-append')
% tau
% error('1')



% taylor series
uppexact = zeros(N+2,1);
upptilde = zeros(N+2,1);
ph = zeros(N+2,1);
xijhat = zeros(N+2,1);
obj.moments;
for i = 2:N+1
   for k = 0:6
   if(mod(k,2) == 0) 
    uppexact(i) = uppexact(i) - (-1)^ceil(k/2) *(2*pi)^(k+2)*sin(2*pi*x(i))*obj.moments(i,k+1)/factorial(k);  
   else
    uppexact(i) = uppexact(i) - (-1)^ceil(k/2) *(2*pi)^(k+2)*cos(2*pi*x(i))*obj.moments(i,k+1)/factorial(k); 
   end
   end
end


% for i = 2:N+1
%     for j = 2:N+1
%         for k = 0:6
%             if(mod(k,2) == 0) 
%                 for l = 0:k
%                     xijhat(k+1) = xijhat(k+1) + nchoosek(k,l)*(x(j)-x(i))^l*obj.moments(j,k-l+1);
% %                                       nchoosek(k,ii-1)*(x1-xi)^(ii-1)*obj.moments(cv1,k-ii+2);  
%                 end
%                 ph(j) = ph(j) - (-1)^ceil(k/2) *(2*pi)^(k+2)*sin(2*pi*x(i))*xijhat(k+1)/factorial(k);               
%             else
%                 for l = 0:k
%                     xijhat(k+1) = xijhat(k+1) + nchoosek(k,l)*(x(j)-x(i))^l*obj.moments(j,k-l+1);
%                     
%                 end
%                 ph(j) = ph(j) - (-1)^ceil(k/2) *(2*pi)^(k+2)*cos(2*pi*x(i))*xijhat(k+1)/factorial(k);               
% 
%             end
%         end
%         upptilde(i) = upptilde(i)+ph(j)*J(i,j);
%     end
% end
% for i = 2:N+1
%     for j = i-1:i+1
%         for k = 0:6
%             if(mod(k,2) == 0) 
%                 for l = 0:k
%                     xijhat(k+1) = xijhat(k+1) + nchoosek(k,l)*(x(j)-x(i))^l*obj.moments(j,k-l+1);
% %                                       nchoosek(k,ii-1)*(x1-xi)^(ii-1)*obj.moments(cv1,k-ii+2);  
%                 end
%                 ph(j) = ph(j) - (-1)^ceil(k/2) *(2*pi)^(k+2)*sin(2*pi*x(i))*xijhat(k+1)/factorial(k);               
%             else
%                 for l = 0:k
%                     xijhat(k+1) = xijhat(k+1) + nchoosek(k,l)*(x(j)-x(i))^l*obj.moments(j,k-l+1);
%                     
%                 end
%                 ph(j) = ph(j) - (-1)^ceil(k/2) *(2*pi)^(k+2)*cos(2*pi*x(i))*xijhat(k+1)/factorial(k);               
% 
%             end
%         end
% %         upptilde(i) = upptilde(i)+ph(j)*J(i,j);
% %         ph(j) = ph(j)+
%     end
%      upptilde(i) = J(i,2:N+1)*ph(2:N+1);
% end

ph;
% obj.moments
upptilde = 0*upptilde;
ph = 0*ph;

for i = 2:N+1
    for j = 2:N+1
        ph(j)=(1/h(j))*(-1/(2*pi))*(cos(2*pi*(x(j)+h(j)/2))-cos(2*pi*(x(j)-h(j)/2)));

    end
    upptilde(i) =  J(i,2:N+1)*ph(2:N+1);
    ph = 0*ph;
end


for j = 2:N+2
   pp(j) = (1/h(j))*(-1/(2*pi))*(cos(2*pi*(x(j)+h(j)/2))-cos(2*pi*(x(j)-h(j)/2))); 
end
[pp' ue];
% max(abs(upptilde-uppexact-tau))
% error('1')
J;

figure
plot(x,uppexact,'*',x,f,'o',x,upptilde,'x')
figure
subplot(223)
plot(x,-(uppexact-upptilde),'^',x,tau,'-*')
ylabel('exact T.E.')
subplot(224)
Y = fft(tau,N)/N;
YY = [Y ifft(Y)];
f = (N)/2*linspace(0,1,N/2+1);

stem(f,2*abs(Y(1:N/2+1))/N) 

% figure
% load('u4.mat')
% tau2 = tau;
% [Z] = obj.unstructuredrecon(u420,p,'solution')
%  [tau]=obj.computefluxintegral(Z,'solution');
%  plot(x,tau2,'+',x,tau,'*')
%  legend('\tau_2','2nd order flux integral of u_4 ')
%  mean(abs(tau-tau2))
%  error('1')

% error('1')
% upptilde
% %   errerr2= NaN;
% %     cverr2 = NaN;
% %     exacterr = NaN;
% %     ee = NaN;
% %     
% %     x10b = x;
% %     h10b = h;
% %     u10b = u;
% %     tau10b = tau;
%     save('te.mat','x10b','h10b','u10b','tau10b','-append')
% 
%     for j = 1:N+2
%         for k = 1:N+2
%           if(abs(J(j,k)) < 1e-4)
%              J(j,k) = 0; 
%           end
%         end
%     end
%     
% errerr2 = nnz(J);
% return 
% %  s = sum(J,2)
% %  t = sum(J',2)
%  for j = 3:N-3:N
%         for k = 1:N+2
%           if(abs(J(j,k)) > 1e-4)
%              J(j,k)=J(j,k)-s(j)/5; 
%           end
%         end
%  end
%     J
% % z = sum(J,2)
% % errerr2 = abs(z(3))+abs(z(end-2));
% % [d,c] = eig(J(2:N+1,2:N+1));
% % c = diag(c);
% % c
% % plot(real(c),imag(c),'*')
% % errerr2 = norm(J-J')
% % return;
% figure
% subplot(221);
% plot(x(2:N+1),ue(2:N+1)-u(2:N+1),'*')
% u = ue(2:N+1)-u(2:N+1);
% Y = fft(u,N)/N;
% f = N/2*linspace(0,1,N/2+1);
% 
% ylabel('exact D.E.')
% subplot(222);
% stem(f,2*abs(Y(1:N/2+1))/N) 
% Y;
% subplot(224)
% Y(1:N/2+1) = Y(1:N/2+1).*-(2*pi*f').^2;
% size(Y(N/2+1:end))
% size(f(end:-1:2))
% 
% Y(N/2+2:end) = Y(N/2+2:end).*-(2*pi*f(end-1:-1:2)').^2;
% f(end:-1:2)
% stem(f,2*abs(Y(1:N/2+1))/N) 
% YY = [ifft(Y) Y]
% uu = real(ifft(Y))*N;
% subplot(223)
% plot(x(2:N+1),uu,'*-')
% ylabel('exact D.E. *(2*pi*f)^2 in frequecy domain')

% u420=u;
% save('u4.mat','-append','u420')
% error('1')

% save('test.mat','tau')
load('test.mat')
% error('1')
    if(q> 0 && r>0)

        fprintf('\n\n\nError Equation:\n')
        obj.computerespseudo();
  
        [Zr] = obj.unstructuredrecon(u,r,'residual');
        figure
        obj.reconplot(Zr,'residual');
  
%     [left,right] = computeflux(Zr,h,N,r,physics,'residual',obj);
%     Rend= (right-left)./h-f;
        Rend = obj.computefluxintegral(Zr,'residual');
 
        
%         figure
%         plot(x,tau,'o',x,-Rend,'*')
%         sum(abs(tau(2:N+1)+Rend(2:N+1)))/N
%         error('1')
        
%    error('2')
 % clearvars -except x obj 
        obj.computeerrorpseudo();

        Zu = obj.unstructuredrecon(obj.convSoln,obj.qOrder,'error');
        obj.convSolnRecon = Zu;
% test using TE for p~=q 
%   load('tauB.mat')
%  [Zq] = obj.unstructuredrecon(u,q,'error');%ue,x,h,N,NaN,NaN,p);
% Rq = obj.computefluxintegral(Zq,'residual');

% [tau4 tau4-Rq -Rend]
% error('1')

%%%%error equation
        uRL = 0;
        uRR = 0;
        for k = 1:p
           uRL = uRL + Z(k,2)*(-h(2)/2)^(k-1);
           uRR = uRR + Z(k,N+1)*(h(N+1)/2)^(k-1);
        end

        if(obj.bcLeftType == 'D')
            obj.bcLeftVal = 0;uRL;(Z(1,2)+Z(2,2)*-h(2)/2); 0;
        end
        if(obj.bcRightType == 'D')
            obj.bcRightVal = 0;uRR;(Z(1,N+1)+Z(2,N+1)*h(N+1)/2);0;
        end

[obj.bcLeftVal obj.bcRightVal];

        exacterr = ue-u;

%  obj.computeerrorpseudo();
        [Z] = obj.unstructuredrecon(exacterr,q,'error');%ue,x,h,N,NaN,NaN,p);

% Z
% figure
% obj.reconplot(Z,'error')
% hold on
% plot(x,exacterr)
% error('1')

% save('test.mat','Rend','-append')
% error('1')
% RR = Rend;
load('test.mat')
% figure
% plot(x(2:N+1),RR(2:N+1),x(2:N+1),Rend(2:N+1))
% error('1')

        f = tau-Rend;
        ff = tau+randn(N+2,1)/N^6;
        ep = rand(N+2,1)/N^6;
         mean(abs(f(2:N+1)-tau(2:N+1)));

         figure
         subplot(1,2,1)
          plot(x,f,'*',x,tau,'o')
          subplot(1,2,2)
          plot(x,tau-f,'x')
          
          f24 = f;
%           xu = x;
%           tauu = tau;
%           save('unstructuredteperiodic.mat','x','tau','f24','-append') 
% error('1')

% % f-fexact test
% p = obj.pOrder;
% obj.pOrder = obj.qOrder;
% obj.computeprimalpseudo();
% obj.computehigherpseudo();
% [ZZ] = obj.unstructuredrecon(obj.exactSolution,obj.pOrder,'solution');
% tauq = obj.computefluxintegral(ZZ,'solution')
% obj.pOrder = p;
%     obj.computeprimalpseudo();
%     obj.computehigherpseudo();
% tau;
% 
% r = obj.rOrder;
% obj.rOrder = obj.qOrder;
% obj.computerespseudo();
% [ZZZ] = obj.unstructuredrecon(u,obj.rOrder,'residual');
% Nr = obj.computefluxintegral(ZZZ,'residual');
% obj.rOrder = r;
% obj.computerespseudo();
% 
% fexact = tauq-Nr;
% % error('1')
% % % % mean(abs(f(2:N+1)-fexact(2:N+1)))
% % % % mean(abs(f(2:N+1)))
% % % % mean(abs(fexact(2:N+1)))
% % % % error('1')

%


% 
% %
% %           error('1')
% F = f;
% % save('tau.mat','x','tau','F') 
% % mean((tau))
% % mean(f)
% % error('1')
% taunew = tefft(x,tau,F);
% % load('tefilt.mat')
% f = taunew;
% 
% % load('tauN.mat')
% % f = tau3264;
% 
% % % load('test.mat')
% % % f = tau;
        obj.errorSource = f;%tau2-Rend;%tau2-Rend;f;%tau6-Rq;%f;%tau;
        
%        [tau2-Rend 2*f]
%        [ tau2 Rend]
%         error('1')
%         z=[tau-ff];
%         [tau ff z]
% 
%         z1 = sum(abs(z(2:N+1)))/N;
%         z2 =sqrt(sum((z(2:N+1)).^2)/N);
%         zinf = max(abs(z(2:N+1)));
%         [z1 z2 zinf]
%         plot(x,tau,x,f)
% 
%         error('1')
 
%    [f -FIr+s]
%    error('1')
 
 
% obj.errorRM
% error('1')

%  [tauE]=reconfluxsoln(Z,f,h,N,q,physics,tlim,obj)
%  error('1')
        [tauE]= obj.computefluxintegral(Z,'error');
        tauE1 = sum(abs(tauE(2:N+1)))/N;
        tauE2 = sqrt(sum((tauE(2:N+1).^2)/N));
        tauEinf = max(abs(tauE(2:N+1)));
%   error('1')

        
        R = ones(N+2,1);
        t=0;
 
        e = exacterr;
        ee = ones(N+2,1);
        ee(1)=NaN;
        ee(N+2) = NaN;
        count = 0;
    
        total_time = 0;
  
        while(max(abs(R)) > 1e-11 )
            Je = obj.computefluxjacobian(e,'error');%,x,h,N,p);
         
%      Jue = obj.computefluxjacobian(u+e,'error');%,x,h,N,p);
%      Ju = obj.computefluxjacobian(u,'error');%,x,h,N,p);
%      Je = Jue-Ju;
            count = count +1;
%      if(count < 50)
%         dt = kk*(40/N)^2; 
        
            Rratio =norm(Rold(2:N+1),2)/norm(R(2:N+1),2); 

            dt = dtold*c2*Rratio;
            K = Je(2:N+1,2:N+1)+eye(N)/dt;
            [Z] = obj.unstructuredrecon(e,q,'error');%u,x,h,N,NaN,NaN,p);

%  [er]=reconplot(x,h,N,p,Z);
            Rold = R;
            [R]=obj.computefluxintegral(Z,'error');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)
 
%             [K\-R(2:N+1) pinv(K)*-R(2:N+1) (K'*K)\(K'*-R(2:N+1)) R(2:N+1)];

%             del = (K'*K)\(K'*-R(2:N+1));%pinv(K)*-R(2:N+1);%K\-R(2:N+1);
% cond(K)
% eig(K)
% R(2:N+1)
%             del = K\-R(2:N+1);
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

%         [x exacterr ee]
max(abs(ue-(u+ee)));
    

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

% w(2:N+1) = w;
% w(1) = NaN;
% w(N+2) = NaN;
% cverr2 = sqrt(sum((ue(2:N+1)-u(2:N+1)).^2)/N);

        ee = exacterr - w;

        errerr = exacterr-ee;
        errerr1 = sum(abs(errerr(2:N+1)))/N;
        errerr2 = sqrt(sum((errerr(2:N+1)).^2)/N) ;
        errerrinf = max(abs(errerr(2:N+1)));

        figure
        plot(x,ee,'*',x,exacterr,'o')

        fprintf('\nError T.E.: [%e\t %e\t %e]\n',tauE1, tauE2, tauEinf);
        fprintf('Error D.E.: [%e\t %e\t %e]\n',errerr1, errerr2, errerrinf);
  
%         ecor = exacterr- (u + ee);
ecor = ue-(u+ee);
% u+ee
%         sqrt(sum((ecor(2:N+1)).^2)/N) 
        
%         Je(2:N+1,2:N+1)
% ee
% J
% Je
    else
%         u
        errerr2 = NaN;
        exacterr = NaN;
        ee = NaN;
    
    end

% tauE

% tau

% J(2:N+1,2:N+1)*exacterr(2:N+1)
% tau(2:N+1)
%  J(2:N+1,2:N+1)*ue(2:N+1)
%  obj.source(2:N+1)
% Je(2:N+1,2:N+1)*w(2:N+1)
% tauE(2:N+1)
% error('1')


% ee = ee- mean(ee(2:N+1))

% % 
% % x = x(2:end-1);
% % u = ee(2:end-1);
% % L = length(x);
% % NFFT = 2^nextpow2(L);
% % Fs = 1*(length(x));
% % Y = fft(u,NFFT)/L;
% % f = Fs/2*linspace(0,1,NFFT/2+1);
% % figure
% % subplot(223)
% % f = reshape(f,size(Y(1:NFFT/2+1)));
% % stem(f,2*abs(Y(1:NFFT/2+1))) 
% % M = max(abs(2*abs(Y(1:NFFT/2+1))));
% % ylim([0 M])
% % xlabel('fft of ee ')
% % subplot(224)
% % size(f);
% % size(Y(1:NFFT/2+1));
% % f = reshape(f,size(Y(1:NFFT/2+1)));
% % stem(f,2*abs(Y(1:NFFT/2+1)).*(2*pi*f).^2) 
% % MM = max(2*abs(Y(1:NFFT/2+1)).*(2*pi*f).^2);
% % ylim([0 3*MM])
% % xlabel('fft of ee * (2\pi f)^2')
% % 
% % 
% % 
% % Y = fft(exacterr(2:end-1),NFFT)/L;
% % subplot(221)
% % stem(f,2*abs(Y(1:NFFT/2+1)))
% % ylim([0 M])
% % 
% % subplot(222)
% % stem(f,2*abs(Y(1:NFFT/2+1)).*(2*pi*f).^2)
% % ylim([0 3*MM])
% % % mean(ee(2:N+1))



% Y


% figure
% subplot(221);
% plot(x,exacterr(2:N+1),'*')
% u = exacterr(2:end-1);
% Y = fft(u,NFFT)/L;
% f = Fs/2*linspace(0,1,NFFT/2+1);
% size(f)
% ylabel('exact D.E.')
% subplot(222);
% stem(f,2*abs(Y(1:NFFT/2+1))/N) 
% Y;
% subplot(224)
% Y(1:NFFT/2+1) = Y(1:NFFT/2+1).*-(2*pi*f').^2;
% size(Y(NFFT/2+1:end))
% size(f(end:-1:2))
% 
% Y(NFFT/2+2:end) = Y(NFFT/2+2:end).*-(2*pi*f(end-1:-1:2)').^2;
% f(end:-1:2)
% stem(f,2*abs(Y(1:NFFT/2+1))/N) 
% YY = [ifft(Y) Y]
% uu = real(ifft(Y))*N;
% subplot(223)
% plot(x,uu,'*-')
% ylabel('exact D.E. *(2*pi*f)^2 in frequecy domain')
end

