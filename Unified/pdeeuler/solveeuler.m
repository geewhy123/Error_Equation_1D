function  [errerr2,x,cverr2,exacterr,ee,te  ] = solveeuler( obj )
%SOLVEEULER Summary of this function goes here
%   Detailed explanation goes here

%errordriver(10,2,0,0,0,'D',[0.904828205821313 0.523900072935957 0.869359219563748],'D',[0.18462798898162 1.854354424142687 0.093932645732845],10,7,'EulerQ','SS');
gam = 1.4;
initializeeuler(obj);
N = obj.nCells
v = zeros(N+2,3);
U = zeros(N+2,3);
x = obj.cellCentroids;
V = obj.initialSolution

% error('1')
p = obj.pOrder;
figure
plot(x,V(:,1),x,V(:,2),x,V(:,3))
%reconstruct
obj.computeprimalpseudo();
Z = obj.unstructuredrecon(V,p,'solution');

figure
obj.reconplot(Z(1:p,:),'solution')
figure
obj.reconplot(Z(p+1:2*p,:),'solution')
figure
obj.reconplot(Z(2*p+1:3*p,:),'solution')
%

% Z
% error('1')

% flux at boundaries
[phi1,phi2,phi3]=computeeulerfluxintegral(obj,Z,'solution');
% error('2')
h = obj.cellWidths(2);
k = 0.00005;
% while(norm([phi1 phi2 phi3]) > 0.1)

% updateboundaryghost(U);

[phi1,phi2,phi3]
% error('1')

   entropy = log(V(:,3)./(V(:,1).^gam));
   figure
   plot(x,entropy,'o')
error('1')


for j = 1:5000

 for i = 2:N+1
 [U(i,1) U(i,2) U(i,3)] = toconservedvars(V(i,1),V(i,2),V(i,3));
 end

   unew = U-k*[phi1 phi2 phi3];
   
  
 for i = 2:N+1
 [V(i,1) V(i,2) V(i,3)] = toprimitivevars(unew(i,1),unew(i,2),unew(i,3));
 end
% % %  updateboundaryghost(unew,h,k,N);
    Z = obj.unstructuredrecon(V,p,'solution');
   

    
% %      if(j==100)
% %    [phi1 phi2 phi3]
% %    rho = Z(1,2)-0.05*Z(2,2); 
% %    u = Z(3,2)-0.05*Z(4,2) ;
% %    P = Z(5,2)-0.05*Z(6,2);
% %    [rho u P]
% %    gam  = 1.4;
% %    
% %    M = u/sqrt(gam*P/rho);
% %    P0 = P*(1+((gam-1)/2)*M^2)^(gam/(gam-1));
% %    T0 = (P/rho)*(1+((gam-1)/2)*M^2);
% %    [T0 P0]
% %    Pb = Z(5,N+1)+0.05*Z(6,N+1)
% %       figure
% % obj.reconplot(Z(1:obj.pOrder,:),'solution')
% % 
% % obj.reconplot(Z(obj.pOrder+1:2*obj.pOrder,:),'solution')
% % 
% % obj.reconplot(Z(2*obj.pOrder+1:3*obj.pOrder,:),'solution')
% %     error('2')
% %    end

    
    
%    figure
% obj.reconplot(Z(1:p,:),'solution')
% 
% obj.reconplot(Z(p+1:2*p,:),'solution')
% 
% obj.reconplot(Z(2*p+1:3*p,:),'solution')
% error('2')
   
   
   [phi1,phi2,phi3]=computeeulerfluxintegral(obj,Z,'solution');
  
   

end
   entropy = log(V(:,3)./(V(:,1).^gam));
   figure
   plot(x,entropy,'o')


figure
% plot(x,unew(:,1),x,unew(:,2),x,unew(:,3))
plot(x,V(:,1),x,V(:,2),x,V(:,3))
[phi1 phi2 phi3]

[mean(phi1(2:N+1)) mean(phi2(2:N+1)) mean(phi3(2:N+1))]


errerr2 = NaN;
cverr2 = NaN;
ee = NaN;
exacterr = NaN;
te = NaN;
end

