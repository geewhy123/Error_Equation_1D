function [ FI ] = computefluxintegral( obj,Z,eqn )
if(strcmp(obj.physics,'Poisson')==1)
    FI=computepoissonfluxintegral(obj,Z,eqn);
elseif(strcmp(obj.physics,'Advection')==1)
    FI=computeadvectionfluxintegral(obj,Z,eqn);
elseif(strcmp(obj.physics,'BurgersMod')==1)
    FI=computeburgersmodfluxintegral(obj,Z,eqn);
else
    assert(0)
end

end

function [ FI ] = computepoissonfluxintegral( obj,Z,eqn )
%COMPUTEFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here
x = obj.cellCentroids;
h = obj.cellWidths;
N = obj.nCells;

if(strcmp(eqn,'solution')==1)
p = obj.pOrder;
elseif(strcmp(eqn,'error')==1)
    p = obj.qOrder;
elseif(strcmp(eqn,'residual')==1)
    p = obj.rOrder;
else
   assert(0); 
end


% Z=unstructuredrecon(u,x,h,N,NaN,NaN,p);

Fr = zeros(N+2,1);
Fl = zeros(N+2,1);
FrAve = zeros(N+2,1);
FlAve = zeros(N+2,1);
F = zeros(N+2,1);
% reconplot(x,h,N,p,Z)

jump = zeros(N+2,1);
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


%       jump = zeros(N+2,1);
      
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

  
  
  %temp
%     for i = 3:N+1
%         ul1 = Z(1,i)+Z(2,i)*(-h(i)/2)+Z(3,i)*(-h(i)/2)^2;
%         ul2 = Z(1,i-1) + Z(2,i-1)*(h(i-1)/2)+Z(3,i-1)*(h(i-1)/2)^2;
%       jump(i) = (.2/((h(i)+h(i-1))/2))*(ul1-ul2) ;
%     
%       
%    end
%     if(obj.bcLeftType=='P' && obj.bcRightType == 'P')
%         i=2;
%         ul1 = Z(1,i)+Z(2,i)*(-h(i)/2);
%         ul2 = Z(1,N+1) + Z(2,N+1)*(h(i-1)/2);
%       jump(i) = (.2/((h(i)+h(i-1))/2))*(ul1-ul2) ;
%     end

  %temp
  
  
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
if(strcmp(eqn,'solution')==1 || strcmp(eqn,'residual')==1)
 FI = (FrAve-FlAve)./h-obj.source;
elseif(strcmp(eqn,'error')==1)
%     size(FrAve)
%     size(obj.errorSource)
  FI = (FrAve-FlAve)./h-obj.errorSource;


end
 
 


% 
%  FI
% error('1')
%  FrAve(3)
% FlAve(3)
% FI(3)
% error('1')

 
end

function [ FI ] = computeadvectionfluxintegral( obj,Z,eqn )
%COMPUTEFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here

% error('1')

if(strcmp(eqn,'solution')==1)
    p = obj.pOrder;
elseif(strcmp(eqn,'error')==1)
    p = obj.qOrder;
elseif(strcmp(eqn,'residual')==1)
    p = obj.rOrder;
else
   assert(0); 
end

x = obj.cellCentroids;
h = obj.cellWidths;
N = obj.nCells;
% p = obj.pOrder;



Fr = zeros(N+2,1);
Fl = zeros(N+2,1);
FrAve = zeros(N+2,1);
FlAve = zeros(N+2,1);


if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')

for i=2:N


    for k = 1:p
   Fr(i) = Fr(i) + Z(k,i+1)*(-h(i+1)/2)^(k-1); 
%    Fl(i) = Fl(i) + k*Z(k+1,i)*(-h(i)/2)^(k-1);
   
    end

end
for k = 1:p
   Fr(N+1) = Fr(N+1) + Z(k,2)*(-h(2)/2)^(k-1); 
   
  
end

    Fr(1) =Fr(N+1) ;
FrAve(2:N+1) = Fr(2:N+1);
FlAve(2:N+1) = Fr(1:N);




elseif(obj.bcLeftType == 'F' && obj.bcRightType == 'D')
for i=2:N

    for k = 1:p
   Fr(i) = Fr(i) + Z(k,i+1)*(-h(i+1)/2)^(k-1); 
%    Fl(i) = Fl(i) + k*Z(k+1,i)*(-h(i)/2)^(k-1);
   
    end

end
for k = 1:p
   Fr(N+1) = Fr(N+1) + Z(k,2)*(-h(2)/2)^(k-1);   
end

    Fr(1) =Fr(N+1) ;
FrAve(2:N+1) = Fr(2:N+1);
FlAve(2:N+1) = Fr(1:N);

% FlAve(2) = 0;
% for k = 1:p
%     FlAve(2) = FlAve(2) + Z(k,2)*(-h(2)/2)^(k-1);   
% end
FrAve(N+1) = obj.bcRightVal;

% Z
% [FlAve FrAve (FrAve-FlAve)./h]

% abs(mean(Fr(2:N+1)))
%  if(abs(mean(Fr(2:N+1))-1) > 1e-12)
% error('1')
%  end

else
   assert(0) 
end


% plot(x,FrAve,x,FlAve)

% % % %  FI = (FrAve-FlAve)./h-obj.source;

 
 if(strcmp(eqn,'solution')==1 || strcmp(eqn,'residual')==1)
 FI = (FrAve-FlAve)./h-obj.source;

 
%  if(strcmp(eqn,'residual')==1)
%      new = [FrAve FlAve]
%  FI
%  error('1')
%  end
 
 
elseif(strcmp(eqn,'error')==1)
%     size(FrAve)
%      size(obj.errorSource)
%      error('1')
  FI = (FrAve-FlAve)./h-obj.errorSource;


end

% 
%  FI
% error('1')
%  FrAve(3)
% FlAve(3)
% FI(3)
% error('1')

 
end



function [ FI ] = computeburgersmodfluxintegral( obj,Z,eqn )
%COMPUTEFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here

% error('1')

if(strcmp(eqn,'solution')==1)
    p = obj.pOrder;
elseif(strcmp(eqn,'error')==1)
    p = obj.qOrder;
elseif(strcmp(eqn,'residual')==1)
    p = obj.rOrder;
else
   assert(0); 
end

x = obj.cellCentroids;
h = obj.cellWidths;
N = obj.nCells;
% p = obj.pOrder;



Fr = zeros(N+2,1);
Fl = zeros(N+2,1);
FrAve = zeros(N+2,1);
FlAve = zeros(N+2,1);


if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')

for i=2:N


    for k = 1:p
   Fr(i) = Fr(i) + Z(k,i+1)*(-h(i+1)/2)^(k-1); 
%    Fl(i) = Fl(i) + k*Z(k+1,i)*(-h(i)/2)^(k-1);
   
    end

end
for k = 1:p
   Fr(N+1) = Fr(N+1) + Z(k,2)*(-h(2)/2)^(k-1); 
   
  
end

    Fr(1) =Fr(N+1) ;
FrAve(2:N+1) = Fr(2:N+1);
FlAve(2:N+1) = Fr(1:N);




elseif(obj.bcLeftType == 'F' && obj.bcRightType == 'D')
 Z
%  obj.reconplot(Z)
 error('1')
    for i=2:N

    for k = 1:p
        
        ur = Z(k,i+1)*(-h(i+1)/2)^(k-1); 
        
%         assert(ur-1 < 0);
   Fr(i) = Fr(i) + ur^2/2-ur; 
%    Fl(i) = Fl(i) + k*Z(k+1,i)*(-h(i)/2)^(k-1);
   
    end

end
for k = 1:p
   ur = Z(k,2)*(-h(2)/2)^(k-1);
%    assert(ur-1 < 0);
   Fr(N+1) = Fr(N+1) + ur^2/2-ur;   
end

    Fr(1) =Fr(N+1) ;
FrAve(2:N+1) = Fr(2:N+1);
FlAve(2:N+1) = Fr(1:N);


ul = obj.bcRightVal;
FrAve(N+1) = ul^2/2-ul;


elseif(obj.bcLeftType == 'D' && obj.bcRightType == 'D')
%      Z
%   obj.reconplot(Z)
%     error('1')
        ur = zeros(N+2,1);
        upr = zeros(N+2,1);
        ul = zeros(N+2,1);
        upl = zeros(N+2,1);

jump = zeros(N+2,1);
for i=2:N+1

    for k = 1:p
        ur(i) = ur(i)+Z(k,i)*(h(i)/2)^(k-1);
        ul(i) = ul(i)+Z(k,i)*(-h(i)/2)^(k-1);
    end
    for k = 1:p-1
        upr(i) = upr(i)+k*Z(k+1,i)*(h(i)/2)^(k-1);
        upl(i) = upl(i)+k*Z(k+1,i)*(-h(i)/2)^(k-1);
    end

%         assert(ur(i)-1 < 0);
   Fr(i) = -1*(ur(i)^2/2-ur(i)-upr(i)); 
   Fl(i) = -1*(ul(i)^2/2-ul(i)-upl(i));
   
   
     if(p==2)
        jump(i) = (.2/((h(i)+h(i-1))/2))*(ul(i)-ur(i-1)) ;
      end


   
end
jump(2) = 0;
for i=2:N+1
   
%         ul(i) = Z(k,i)*(-h(i+1)/2)^(k-1);    
%    Fl(i) = Fl(i) + k*Z(k+1,i)*(-h(i)/2)^(k-1);
   



    if i==2
        FrAve(i) = (Fr(i)+Fl(i+1))/2;
        FlAve(i) = Fl(i); 
    
        
    elseif i==N+1
        FrAve(i) = Fr(i);
        FlAve(i) = (Fr(i-1)+Fl(i))/2;
    else
%         jumpr = 0;
%         jumpl = 0;
      
        
        
      FrAve(i) = (Fr(i)+Fl(i+1))/2+jump(i+1);
      FlAve(i) = (Fr(i-1)+Fl(i))/2+jump(i);
      

    end



end

% xx = 0:0.1:1;
% ff = -0.5*sech(xx/2).^2-0.5*(1-tanh(xx/2).^2)+(1-tanh(xx/2));
%     [xx' ff']

%      [ul ur jump]
%    error('1')
% [ul ur upl upr Fl Fr FlAve FrAve]
% error('1')
    
% for k = 1:p
%    ur = Z(k,2)*(-h(2)/2)^(k-1);
% %    assert(ur-1 < 0);
%    Fr(N+1) = Fr(N+1) + ur^2/2-ur;
%    
% 
%    
% end
% 
%     Fr(1) =Fr(N+1) ;
% FrAve(2:N+1) = Fr(2:N+1);
% FlAve(2:N+1) = Fr(1:N);
% 
% 
% ur = obj.bcRightVal;
% FrAve(N+1) = ur^2/2-ur;
% 
% 
% ul = obj.bcLeftVal;
% FlAve(2) = ul^2/2-ul;
% [Fr FlAve FrAve]
% error('1')


else
   assert(0) 
end


 
 if(strcmp(eqn,'solution')==1 || strcmp(eqn,'residual')==1)

% [FlAve FrAve]
% error('2')
     
     FI = (FrAve-FlAve)./h-obj.source;

 
 
elseif(strcmp(eqn,'error')==1)

  FI = (FrAve-FlAve)./h-obj.errorSource;


end

 
end