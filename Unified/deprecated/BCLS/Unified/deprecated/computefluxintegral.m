function [ FI ] = computefluxintegral( u,x,h,N,p )
%COMPUTEFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here
error('1')
Z=unstructuredrecon(u,x,h,N,NaN,NaN,p);

Fr = zeros(N+2,1);
Fl = zeros(N+2,1);
FrAve = zeros(N+2,1);
FlAve = zeros(N+2,1);
F = zeros(N+2,1);
reconplot(x,h,N,p,Z)

for i=2:N+1
    
    for k = 1:p-1
   Fr(i) = Fr(i) + k*Z(k+1,i)*(+h(i)/2)^(k-1); 
   Fl(i) = Fl(i) + k*Z(k+1,i)*(-h(i)/2)^(k-1);
   
    end


    
end


for i=2:N+1
    
    if i==2
        FrAve(i) = (Fr(i)+Fl(i+1))/2;
        FlAve(i) = (Fr(N+1)+Fl(i))/2;
        
    elseif i==N+1
        FrAve(i) = (Fr(i)+Fl(2))/2;
        FlAve(i) = (Fr(i-1)+Fl(i))/2;
    else
      FrAve(i) = (Fr(i)+Fl(i+1))/2;
      FlAve(i) = (Fr(i-1)+Fl(i))/2;
    end
    
    
end


% figure
% plot(x,Fr,x,Fl);
% FrAve(N+1) = Fr(N+1);
% FlAve(2) = Fl(2);

Fr;
Fl;

% plot(x,FrAve,x,FlAve)
 FI = (FrAve-FlAve)./h;

end

