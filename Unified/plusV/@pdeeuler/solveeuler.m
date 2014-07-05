function  [errerr2,x,cverr2,exacterr,ee,te  ] = solveeuler( obj )
%SOLVEEULER Summary of this function goes here
%   Detailed explanation goes here


obj.T0 = obj.bcLeftVal(2);
obj.P0 = obj.bcLeftVal(1);
obj.Pb = obj.bcRightVal(1);
obj.primalT0 = obj.T0;
obj.primalP0 = obj.P0;
obj.primalPb = obj.Pb;
obj.areatype = obj.bcRightVal(2);

%errordriver(10,2,0,0,0,'D',[0.904828205821313 0.523900072935957 0.869359219563748],'D',[0.18462798898162 1.854354424142687 0.093932645732845],10,7,'EulerQ','SS');

obj.initializeeuler();
N = obj.nCells;
v = zeros(N+2,3);
U = zeros(N+2,3);
x = obj.cellCentroids;
V = obj.initialSolution;

p = obj.pOrder;
q = obj.qOrder;
r = obj.rOrder
figure
plot(x,V(:,1),x,V(:,2),x,V(:,3))
%reconstruct
obj.computeprimalpseudo();



%%%
% obj
% 
if(strcmp(obj.goal,'SS')==1 && 0)
    [errerr2,x,cverr2,exacterr,ee,te  ] =obj.solvebyeulerjacobian();
    U = obj.convSoln;
    for i = 2:N+1
        [ V(i,1),V(i,2),V(i,3)] = toprimitivevars(U(i,1),U(i,2),U(i,3)); 
    end
else

%%%explicit time stepping

    Z = obj.unstructuredrecon(V,p,'solution');

%     figure
%     obj.reconplot(Z(1:p,:),'solution')
%     figure
%     obj.reconplot(Z(p+1:2*p,:),'solution')
%     figure
%     obj.reconplot(Z(2*p+1:3*p,:),'solution')
%     %




% flux at boundaries
[phi1,phi2,phi3]=obj.computeeulerfluxintegral(Z,'solution');
% error('2')
h = obj.cellWidths;
C = 0.6;
k = C*mean(h(2:N+1));%0.02;
% while(norm([phi1 phi2 phi3]) > 0.1)

% updateboundaryghost(U);

[phi1,phi2,phi3]

for i = 2:N+1
[U(i,1),U(i,2),U(i,3)]=toconservedvars(V(i,1),V(i,2),V(i,3));
end

% error('1')

[mean(abs(phi1(2:N+1))) mean(abs(phi2(2:N+1))) mean(abs(phi3(2:N+1)))]
%  error('1')


count = 0;
while(max(max( max(abs(phi1(2:N+1)),abs(phi2(2:N+1))),abs(phi3(2:N+1)))) > 1e-3 )
count = count+1;%floor(100/k)%5000

 for i = 2:N+1
 [U(i,1) U(i,2) U(i,3)] = toconservedvars(V(i,1),V(i,2),V(i,3));
 end

   unew = U-k*[phi1 phi2 phi3];
   
  
 for i = 2:N+1
 [V(i,1) V(i,2) V(i,3)] = toprimitivevars(unew(i,1),unew(i,2),unew(i,3));
 end
% % %  updateboundaryghost(unew,h,k,N);
    Z = obj.unstructuredrecon(V,p,'solution');
      
   [phi1,phi2,phi3]=obj.computeeulerfluxintegral(Z,'solution');
  
   if(mod(count,100)==0)
[mean(abs(phi1(2:N+1))) mean(abs(phi2(2:N+1))) mean(abs(phi3(2:N+1)))]
   end
end
u = unew;

exacterrv = obj.exactSolutionV-V

cverr2 = [sqrt(sum((exacterrv(2:N+1,1)).^2)/N) sqrt(sum((exacterrv(2:N+1,2)).^2)/N) sqrt(sum((exacterrv(2:N+1,3)).^2)/N)]


%
%te
    Z = obj.unstructuredrecon(obj.exactSolutionV,p,'solution');
[tauu1 tauu2 tauu3]=obj.computeeulerfluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,tlim,obj)
  [tauu1 tauu2 tauu3]
te1 = [ sum(abs(tauu1(2:N+1)))/N  sum(abs(tauu2(2:N+1)))/N sum(abs(tauu3(2:N+1)))/N ]
teu = [tauu1 tauu2 tauu3];
 te = teu;



load('test.mat')
fprintf('using saved info') 

 
%residual
 obj.convSoln = u;
 u=obj.convSoln;
obj.convSolutionV = V;
 [Z] = obj.unstructuredrecon(V,p,'solution');
obj.convVreconp = Z;



%  save('test.mat','u','V','Z','te');


obj.computeprimalleftright();



if(q> 0 && r>0)

  obj.computerespseudo();
   
  
  Zr = obj.unstructuredrecon(V,r,'residual');



  
  
[R1, R2, R3] = obj.computeeulerfluxintegral(Zr,'residual');
 
[R1 R2 R3]

norm1R  = [sum(abs(R1(2:N+1)))/N sum(abs(R2(2:N+1)))/N sum(abs(R3(2:N+1)))/N]


obj.computeerrorpseudo();



 Zu = obj.unstructuredrecon(V,obj.qOrder,'error');
 obj.convSolnRecon = Zu;
 
 
 

%%%error equation

if(obj.bcLeftType == 'D')
   obj.T0 = 0; 
   obj.P0 = 0;
end
if(obj.bcRightType == 'D')
    obj.Pb = 0;
end

exacterrv = obj.exactSolutionV-V;
exacterru = obj.exactSolutionU-u

obj.exactSolutionV
UU =zeros(N+2,3);
for z = 2:N+1
[UU(z,1),UU(z,2),UU(z,3)] = toconservedvars(obj.exactSolutionV(z,1),obj.exactSolutionV(z,2),obj.exactSolutionV(z,3)); 
end
UU

obj.exactSolutionU


f = -[R1 R2 R3];
   obj.errorSource = teu;%f;%tau;
   fprintf('using t.e. as source')
 
 



 del = ones(3*N,1);
 R = ones(3*N+2,1);
 t=0;
 
 e = exacterrv;
  ev = exacterrv;

 ee = ones(3*N+2,1);
 ee(1)=NaN;
 ee(3*N+2) = NaN;
 count = 0;
 E = NaN*ones(3*N+2,1);
  
dtold = 0.01;
  c2 = 10;
k = k/2;

Upe =zeros(N+2,3);
Vpe = zeros(N+2,3);
Upenew =zeros(N+2,3);
Vpenew = zeros(N+2,3); 

  Z = obj.unstructuredrecon(e,q,'error');
[phi1,phi2,phi3]=obj.computeeulerfluxintegral(Z,'error');
% V
% obj.convVreconp
% [obj.convVleft obj.convVright]
%  [phi1 phi2 phi3]
% error('1')
max(max( max(abs(phi1(2:N+1)),abs(phi2(2:N+1))),abs(phi3(2:N+1))))
while(max(max( max(abs(phi1(2:N+1)),abs(phi2(2:N+1))),abs(phi3(2:N+1)))) > 1e-13 )
count = count+1;%floor(100/k)%5000


 for i = 2:N+1
 Vpe(i,1) = V(i,1)+e(i,1);
 Vpe(i,2) = V(i,2)+e(i,2);
 Vpe(i,3) = V(i,3)+e(i,3);
 end


 for i = 2:N+1
 [Upe(i,1) Upe(i,2) Upe(i,3)] = toconservedvars(Vpe(i,1),Vpe(i,2),Vpe(i,3));
 end


 E = Upe-u;
 
   enew = E-k*([phi1 phi2 phi3]);
   
  Upenew = u+enew;
   
 for i = 2:N+1
 [Vpenew(i,1) Vpenew(i,2) Vpenew(i,3)] = toprimitivevars(Upenew(i,1),Upenew(i,2),Upenew(i,3));
 end

 for i = 2:N+1
 ev(i,1) = Vpenew(i,1)-V(i,1);
 ev(i,2) = Vpenew(i,2)-V(i,2);
 ev(i,3) = Vpenew(i,3)-V(i,3);
 end

  
% % %  updateboundaryghost(unew,h,k,N);
    Z = obj.unstructuredrecon(ev,q,'error');
    
   [phi1,phi2,phi3]=obj.computeeulerfluxintegral(Z,'error');
   e = ev;


   if(mod(count,100)==0)
[mean(abs(phi1(2:N+1))) mean(abs(phi2(2:N+1))) mean(abs(phi3(2:N+1)))]
enew
   end
end



% while(max(abs(R(2:3*N+1))) > 1e-13  )
% %     Je = obj.computeeulerfluxjacobian(e,'error');%,x,h,N,p);
%      count = count +1;
%         
% %      Rratio =norm(Rold(2:3*N+1),2)/norm(R(2:3*N+1),2); 
% 
%      dt= dtold*c2*Rratio;
% 
%      K = Je(2:3*N+1,2:3*N+1)+eye(3*N)/dt;
% 
% 
%     [Z] = obj.unstructuredrecon(e,q,'error');%u,x,h,N,NaN,NaN,p);
% 
% 
%     Rold = R;
%     [phi1 phi2 phi3]=obj.computeeulerfluxintegral(Z,'error');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)
%     
%     R(2:3:3*N-1) = phi1(2:N+1);
%     R(3:3:3*N) = phi2(2:N+1);
%     R(4:3:3*N+1) = phi3(2:N+1);
%     eu = NaN*ones(N+2,3);
% 
% 
% %%%%%new translation
%     Vpe = e+V;
% 
% 
%     for j = 2:N+1
% %     eu(j,1) = e(j,1);
% %     eu(j,2) = e(j,2);
% %     eu(j,3) = e(j,3);
% % [eu(j,1),eu(j,2),eu(j,3)] = toconservedvars(e(j,1),e(j,2),e(j,3));
%         [Upe(j,1),Upe(j,2),Upe(j,3)] = toconservedvars(Vpe(j,1),Vpe(j,2),Vpe(j,3));
%     end
% 
%     eu = Upe-u;
% % error('1')
%  
% 
% %%%%%new translation
%  
% 
%     del = K\-R(2:3*N+1);
% 
%     [Rratio max(abs(R(2:3*N+1)))]
% 
%     
%     E(2:3:3*N-1) = eu(2:N+1,1);
%     E(3:3:3*N) = eu(2:N+1,2);
%     E(4:3:3*N+1) = eu(2:N+1,3);
% 
% 
%      EE = E(2:3*N+1) + del;%*dt;
%      E = NaN*ones(N+2,1);
%      E(2:3*N+1) = EE;
%      
%      eu(2:N+1,1) = E(2:3:3*N-1) ;
%      eu(2:N+1,2) = E(3:3:3*N) ;
%      eu(2:N+1,3) = E(4:3:3*N+1);
% 
% 
% %%%%%%new trans
% 
%     Upe = eu+u;
% 
% 
%     for j = 2:N+1
% %     e(j,1) = eu(j,1);
% %     e(j,2) = eu(j,2);
% %     e(j,3) = eu(j,3);
% % [e(j,1),e(j,2),e(j,3)] = toprimitivevars(eu(j,1),eu(j,2),eu(j,3));
%         [Vpe(j,1),Vpe(j,2),Vpe(j,3)] = toprimitivevars(Upe(j,1),Upe(j,2),Upe(j,3));
%     end
%      
%     
%     e = Vpe-V;
% % error('1')
% %%%%%%new trans
% 
% 
% %      e
% % exacterr
% 
% % error('1')
%      
%      t = t+dt;
%      dtold = dt;
%  end

ee = e;
ee

w = exacterrv-ee

count

max(abs(w))
% error('1')
figure
plot(x,w)



ee = exacterrv - w;
errerr2 = [sqrt(sum((exacterrv(2:N+1,1)-ee(2:N+1,1)).^2)/N) sqrt(sum((exacterrv(2:N+1,2)-ee(2:N+1,2)).^2)/N)  sqrt(sum((exacterrv(2:N+1,3)-ee(2:N+1,3)).^2)/N) ]
w
figure
subplot(3,1,1)
plot(x,ee(:,1),'*',x,exacterrv(:,1),'o')
subplot(3,1,2)
plot(x,ee(:,2),'*',x,exacterrv(:,2),'o')
subplot(3,1,3)
plot(x,ee(:,3),'*',x,exacterrv(:,3),'o')

exacterr = exacterrv
else
   errerr2 = NaN;
   exacterr = NaN;
   ee = NaN;
    
end




end





end
