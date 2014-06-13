function  [errerr2,x,cverr2,exacterr,ee,te  ] = solveeuler( obj )
%SOLVEEULER Summary of this function goes here
%   Detailed explanation goes here

%errordriver(10,2,0,0,0,'D',[0.904828205821313 0.523900072935957 0.869359219563748],'D',[0.18462798898162 1.854354424142687 0.093932645732845],10,7,'EulerQ','SS');

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



% flux at boundaries
[phi1,phi2,phi3]=computeeulerfluxintegral(obj,Z,'solution')
% error('2')
h = obj.cellWidths(2);
k = 0.001;
% while(norm([phi1 phi2 phi3]) > 0.1)

% updateboundaryghost(U);



for j = 1:10000

 for i = 2:N+1
 [U(i,1) U(i,2) U(i,3)] = toconservedvars(V(i,1),V(i,2),V(i,3));
 end

   unew = U+k*[phi1 phi2 phi3];
   
   [phi1 phi2 phi3]
%     error('2')

 for i = 2:N+1
 [V(i,1) V(i,2) V(i,3)] = toprimitivevars(unew(i,1),unew(i,2),unew(i,3));
 end
 updateboundaryghost(unew,h,k,N);
    Z = obj.unstructuredrecon(V,p,'solution');
   

%    figure
% obj.reconplot(Z(1:p,:),'solution')
% 
% obj.reconplot(Z(p+1:2*p,:),'solution')
% 
% obj.reconplot(Z(2*p+1:3*p,:),'solution')
% error('2')
   
   
   [phi1,phi2,phi3]=computeeulerfluxintegral(obj,Z,'solution');
   
end
figure
plot(x,unew(:,1),x,unew(:,2),x,unew(:,3))
[phi1 phi2 phi3]


errerr2 = NaN;
cverr2 = NaN;
ee = NaN;
exacterr = NaN;
te = NaN;
end

