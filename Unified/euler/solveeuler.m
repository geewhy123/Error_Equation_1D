function  [errerr2,x,cverr2,exacterr,ee  ] = solveeuler( obj )
%SOLVEEULER Summary of this function goes here
%   Detailed explanation goes here

%errordriver(10,2,0,0,0,'D',[0.904828205821313 0.523900072935957 0.869359219563748],'D',[0.18462798898162 1.854354424142687 0.093932645732845],10,7,'EulerQ','SS');

initializeeuler(obj);


x = obj.cellCentroids;
U = obj.initialSolution
p = obj.pOrder;
figure
plot(x,U(:,1),x,U(:,2),x,U(:,3))
%reconstruct
obj.computeprimalpseudo();
Z = obj.unstructuredrecon(U,p,'solution')

obj.reconplot(Z(1:p,:),'solution')
obj.reconplot(Z(p+1:2*p,:),'solution')
obj.reconplot(Z(2*p+1:3*p,:),'solution')
%


errerr2 = NaN;
cverr2 = NaN;
ee = NaN;
exacterr = NaN;
end

