function [ FI ] = computefluxintegral( obj,Z )
%COMPUTEFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here
x = obj.cellCentroids;
h = obj.cellWidths;
N = obj.nCells;
p = obj.pOrder;



% Z=unstructuredrecon(u,x,h,N,NaN,NaN,p);

Fr = zeros(N+2,1);
Fl = zeros(N+2,1);
FrAve = zeros(N+2,1);
FlAve = zeros(N+2,1);
F = zeros(N+2,1);
% reconplot(x,h,N,p,Z)


for i=2:N+1

% if(i==2)
%    Z(2,i)
% %    error('1')
% end


    for k = 1:p-1
   Fr(i) = Fr(i) + k*Z(k+1,i)*(+h(i)/2)^(k-1); 
   Fl(i) = Fl(i) + k*Z(k+1,i)*(-h(i)/2)^(k-1);
   
    end


    
end


      jump = zeros(N+2,1);
      
  if(p==2)
   for i = 3:N+1
       
    
        ul1 = Z(1,i)+Z(2,i)*(-h(i)/2);
        ul2 = Z(1,i-1) + Z(2,i-1)*(h(i-1)/2);
      jump(i) = (.2/((h(i)+h(i-1))/2))*(ul1-ul2) ;
    
      
   end
    if(obj.bcLeftType=='P' && obj.bcRightType == 'P')
        i=2;
        ul1 = Z(1,i)+Z(2,i)*(-h(i)/2);
        ul2 = Z(1,N+1) + Z(2,N+1)*(h(i-1)/2);
      jump(i) = (.2/((h(i)+h(i-1))/2))*(ul1-ul2) ;
    end
  end

%  
%   jump(:) = 0;

for i=2:N+1
    
    if i==2
        FrAve(i) = (Fr(i)+Fl(i+1))/2+jump(i+1);
        if(obj.bcLeftType == 'P')
        FlAve(i) = (Fr(N+1)+Fl(i))/2+jump(2);
        elseif(obj.bcLeftType =='D')
            FlAve(i) = Fl(i);
        else
            assert(0)
        end
        
    elseif i==N+1
        if(obj.bcRightType == 'P')
        FrAve(i) = (Fr(i)+Fl(2))/2+jump(2);
        elseif(obj.bcRightType=='D')
        FrAve(i) = Fr(i);
        else
            assert(0)
        end
        FlAve(i) = (Fr(i-1)+Fl(i))/2+jump(i);
    else
      FrAve(i) = (Fr(i)+Fl(i+1))/2+jump(i+1);
      FlAve(i) = (Fr(i-1)+Fl(i))/2+jump(i);
    end
    
    
end


% figure
% plot(x,Fr,x,Fl);
% FrAve(N+1) = Fr(N+1);
% FlAve(2) = Fl(2);

Fr;
Fl;





% plot(x,FrAve,x,FlAve)

 FI = (FrAve-FlAve)./h-obj.source;


% 
%  FI
% error('1')
%  FrAve(3)
% FlAve(3)
% FI(3)
% error('1')

 
end

