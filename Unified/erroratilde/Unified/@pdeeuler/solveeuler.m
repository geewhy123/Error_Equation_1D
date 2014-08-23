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
% figure
% plot(x,V(:,1),x,V(:,2),x,V(:,3))

%reconstruct
obj.computeprimalpseudo();


%%%
% obj
% 
if(strcmp(obj.goal,'SS')==1)
    [errerr2,x,cverr2,exacterr,ee,te  ] =obj.solvebyeulerjacobian();
    U = obj.convSoln;
    for i = 2:N+1
        [ V(i,1),V(i,2),V(i,3)] = toprimitivevars(U(i,1),U(i,2),U(i,3)); 
    end
else

%%%explicit time stepping

    Z = obj.unstructuredrecon(V,p,'solution');

    figure
    obj.reconplot(Z(1:p,:),'solution')
    figure
    obj.reconplot(Z(p+1:2*p,:),'solution')
    figure
    obj.reconplot(Z(2*p+1:3*p,:),'solution')
    %


% flux at boundaries
[phi1,phi2,phi3]=obj.computeeulerfluxintegral(Z,'solution');
% error('2')
h = obj.cellWidths;
C = 0.2;
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


for j = 1:5000%floor(100/k)%5000

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
  
   
end


[phi1 phi2 phi3]

[mean(abs(phi1(2:N+1))) mean(abs(phi2(2:N+1))) mean(abs(phi3(2:N+1)))]
end





end

