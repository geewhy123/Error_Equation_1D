function [ FI ] = computefluxintegral( obj,Z,eqn )
if(strcmp(obj.physics,'Poisson')==1)
%     FI=computepoissonfluxintegral(obj,Z,eqn);
    h = obj.cellWidths;
    N = obj.nCells;
    F = zeros(N+2,1);
    FI = zeros(N+2,1);
    if(strcmp(eqn,'solution')==1 || strcmp(eqn,'residual')==1)
    f = obj.source;
    elseif(strcmp(eqn,'error')==1)
        f = obj.errorSource;
        
    end
    for i = 1:N+1
        F(i) = computepoissonflux(obj,Z(:,i),Z(:,i+1),eqn,i);
        
        if(i==1)
            if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')
                F(i) = computepoissonflux(obj,Z(:,N+1),Z(:,i+1),eqn,i);
            end
        elseif(i==N+1)
            if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')
                F(i) = computepoissonflux(obj,Z(:,i),Z(:,2),eqn,i);
            end
        end
        
        if(i>1)
        FI(i) = (F(i)-F(i-1))/h(i)-f(i);
        end
    end
    F
    
    
    
elseif(strcmp(obj.physics,'Advection')==1)
    FI=computeadvectionfluxintegral(obj,Z,eqn);
elseif(strcmp(obj.physics,'BurgersMod')==1)
    FI=computeburgersmodfluxintegral(obj,Z,eqn);
elseif(strcmp(obj.physics,'BurgersVisc')==1)
    if(strcmp(eqn,'error')==1)
    FI=computeburgersviscfluxintegralb(obj,Z,eqn);
    else
    FI=computeburgersviscfluxintegral(obj,Z,eqn);
    end
else
    assert(0)
end

end

function [ F ] = computepoissonflux( obj,left,right,eqn,i )
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



if(~isempty(obj.refinecells))
    if(i==1)
        p = obj.hOrder;
    elseif(i==N+1)
        p = obj.hOrder;
    else
        n = length(obj.refinecells);
        if(i <= obj.refinecells(n/2)) 
            pl = obj.hOrder;
            pr = obj.hOrder;
        elseif(i==obj.refinecells(n/2)+1)
            pl = obj.hOrder;
            pr = obj.pOrder;
        elseif(i== obj.refinecells(n/2+1)+1)
            pl = obj.pOrder;
            pr = obj.hOrder;
        elseif(i>= obj.refinecells(n/2+1))
            pl = obj.hOrder;
            pr = obj.hOrder;
        else
             pl = obj.pOrder;
        pr = obj.pOrder;
        end
            
    end
else
    p = obj.pOrder;
    pr = obj.pOrder;
    pl = orj.pOrder;
end       
 
%     case 1
%         p =obj.hOrder;
%         
%     case 2
%         pl = obj.hOrder;
%         pr = obj.hOrder;
%     case 3
%         pl = obj.hOrder;
%         pr = obj.pOrder;
%     case N-1 
%         pl = obj.pOrder;
%         pr = obj.hOrder;
%     case N
%         pl = obj.hOrder;
%         pr = obj.hOrder;
%     case N+1
%         p = obj.hOrder;
%     otherwise 
%         pl = obj.pOrder;
%         pr = obj.pOrder;
% end


if(i==1 && obj.bcLeftType == 'D' && obj.bcRightType == 'D')
   F = 0;
    for k = 1:p-1
     F = F + k*right(k+1)*(-h(i+1)/2)^(k-1);
    end
    return;
elseif(i==N+1 && obj.bcLeftType == 'D' && obj.bcRightType == 'D')
   F = 0;
    for k = 1:p-1
     F = F + k*left(k+1)*(h(i)/2)^(k-1);
    end
    return;
else
    Fl = 0;
    Fr = 0;
    F = 0;
    for k = 1:pr-1
     Fr = Fr + k*right(k+1)*(-h(i+1)/2)^(k-1);
    end
    for k = 1:pl-1
     Fl = Fl + k*left(k+1)*(h(i)/2)^(k-1);
    end
     F = 0.5*(Fr+Fl);
     
     if(pr==2 && pl == 2)
%          if(i==1)
%              if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')
%                 ul1 = right(1)+right(2)*(-h(i+1)/2);
%                 ul2 = left(1)+ left(2)*(h(N+1)/2);
%                 jump = (.2/((h(i)+h(N+1))/2))*(ul1-ul2) ;
%              end
%          elseif(i==N+1)
%                 ul1 = right(1)+right(2)*(-h(i)/2);
%                 ul2 = left(1)+ left(2)*(h()/2);
%                 jump = (.2/((h(i)+h(N+1))/2))*(ul1-ul2) ;
%              if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')
%                  
%              end
%          
%          else
              ul1 = right(1)+right(2)*(-h(i+1)/2);
            ul2 = left(1) + left(2)*(h(i)/2);
                  jump = (.2/((h(i+1)+h(i))/2))*(ul1-ul2) ;
%          end

      F = F+jump;
     end
%     end
    
    return;
    
end
    
    






Fr = zeros(N+2,1);
Fl = zeros(N+2,1);
FrAve = zeros(N+2,1);
FlAve = zeros(N+2,1);
F = zeros(N+2,1);
% reconplot(x,h,N,p,Z)

jump = zeros(N+2,1);
for i=2:N+1

    
    %%%%higher
if((strcmp(eqn,'solution')==1 && obj.hOrder > obj.pOrder) || (strcmp(eqn,'residual')==1 && obj.hOrder > obj.rOrder) || (strcmp(eqn,'error')==1 && obj.hOrder > obj.qOrder) )
if(i==2 || i == 3  || i == N || i == N+1)
    p(i) = obj.hOrder;
else
   p(i) = obj.pOrder; 

end
end
% %%%%higher

    
    
    

    for k = 1:p(i)-1
   Fr(i) = Fr(i) + k*Z(k+1,i)*(+h(i)/2)^(k-1); 
   Fl(i) = Fl(i) + k*Z(k+1,i)*(-h(i)/2)^(k-1);
 
    end


    
end

[Fl Fr]
% error('1')

%       jump = zeros(N+2,1);
      
%   if(p==2)
   for i = 3:N+1
       if(p(i) == 2)
        ul1 = Z(1,i)+Z(2,i)*(-h(i)/2);
        ul2 = Z(1,i-1) + Z(2,i-1)*(h(i-1)/2);
      jump(i) = (.2/((h(i)+h(i-1))/2))*(ul1-ul2) ;
    ul1
    ul2
    error('1')
       end
   end
   
   
    if(obj.bcLeftType=='P' && obj.bcRightType == 'P')
        i=2;
        if(p(i)==2)
        ul1 = Z(1,i)+Z(2,i)*(-h(i)/2);
        ul2 = Z(1,N+1) + Z(2,N+1)*(h(i-1)/2);
      jump(i) = (.2/((h(i)+h(i-1))/2))*(ul1-ul2) ;
        end
    end
%   end

  
  
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
[FlAve FrAve]
% error('1')


% plot(x,FrAve,x,FlAve)
if(strcmp(eqn,'solution')==1 || strcmp(eqn,'residual')==1)
    
    xx(1) = 0;
    for m = 2:N+1
       xx(m) = obj.cellCentroids(m+1)-obj.cellWidths(m+1)/2;
    end
%     xx
%     error('4')
    
% % %     yy = pi*cos(pi*xx);
% % % %     yy = 2*xx-3*xx.^2;
% % %     yy(end+1) = 0;
% % %     obj.reconplot(Z,eqn)
% % %     size(FlAve)
% % %     size(FrAve)
% % %     size(yy')
% % %     [FlAve FrAve yy' ]
% % % FrAve(2:N)-yy(2:N)'
% % %      cverr1 = sum(abs(FrAve(2:N)-yy(2:N)'))/(N-1)
% % %      abs(Fr(4)-yy(4))
% % %      [Fl Fr]

%     error('2')
    
 FI = (FrAve-FlAve)./h-obj.source;
elseif(strcmp(eqn,'error')==1)

  FI = (FrAve-FlAve)./h-obj.errorSource;


end
 
 
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
     
     ur = zeros(N+2,1);
        upr = zeros(N+2,1);
        ul = zeros(N+2,1);
        upl = zeros(N+2,1);
        utilder = zeros(N+2,1);
        utildel = zeros(N+2,1);
        
        nonlinearerror = 1;
         if(nonlinearerror && strcmp(eqn,'error')==1)
%              obj.computerespseudo();%
%              Zu = obj.unstructuredrecon(obj.convSoln,obj.qOrder,'error');
                      Zu = obj.convSolnRecon;%obj.unstructuredrecon(obj.convSoln,obj.qOrder,'error');
                      uorder = obj.qOrder;
%              obj.computeerrorpseudo();%
%                 obj.convSoln
%                 error('1')
 Zu;
% error('1')
         end
        
         
         
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
   
   
   
   if(nonlinearerror && strcmp(eqn,'error')==1)
        for k = 1:uorder
        utilder(i) = utilder(i)+Zu(k,i)*(h(i)/2)^(k-1);
        utildel(i) = utildel(i)+Zu(k,i)*(-h(i)/2)^(k-1);
        end
  
%           if(i==N+1)
%       [utildel.*ul utilder.*ur Fl Fr]
%       [ul ur]
%       error('1')
%    end
      Fr(i) = Fr(i) - ur(i) * utilder(i);
      Fl(i) = Fl(i) - ul(i) * utildel(i);
      
      
%       if(i==12)
%          ur(i) * utilder(i)
%          Fr(i)
%          error('1')
%       end
  
%     end
    
   end
   
   
   
   
     if(p==2)
        jump(i) = (.2/((h(i)+h(i-1))/2))*(ul(i)-ur(i-1)) ;
      end


   
end

% error('1')


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



function [ FI ] = computeburgersviscfluxintegral( obj,Z,eqn )
%COMPUTEFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here

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


Fr = zeros(N+2,1);
Fl = zeros(N+2,1);
FrAve = zeros(N+2,1);
FlAve = zeros(N+2,1);


if(obj.bcLeftType == 'D' && obj.bcRightType == 'D')
     
     ur = zeros(N+2,1);
        upr = zeros(N+2,1);
        ul = zeros(N+2,1);
        upl = zeros(N+2,1);
        utilder = zeros(N+2,1);
        utildel = zeros(N+2,1);
        
        nonlinearerror = 1;
         if(nonlinearerror && strcmp(eqn,'error')==1)
%              obj.computerespseudo();%
%              Zu = obj.unstructuredrecon(obj.convSoln,obj.qOrder,'error');
                      Zu = obj.convSolnRecon;%obj.unstructuredrecon(obj.convSoln,obj.qOrder,'error');
                      uorder = obj.qOrder;
%              obj.computeerrorpseudo();%
%                 obj.convSoln
%                 error('1')
 Zu;
% error('1')
         end
                
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
   Fr(i) = -1*(ur(i)^2/2-upr(i)); %factor the flux function?
   Fl(i) = -1*(ul(i)^2/2-upl(i));
    
   if(nonlinearerror && strcmp(eqn,'error')==1)
        for k = 1:uorder
        utilder(i) = utilder(i)+Zu(k,i)*(h(i)/2)^(k-1);
        utildel(i) = utildel(i)+Zu(k,i)*(-h(i)/2)^(k-1);
        end
  

      Fr(i) = Fr(i) - ur(i) * utilder(i);
      Fl(i) = Fl(i) - ul(i) * utildel(i);
      
   end

     if(p==2)
        jump(i) = (.2/((h(i)+h(i-1))/2))*(ul(i)-ur(i-1)) ;
      end

end

jump(2) = 0;
for i=2:N+1


    if i==2
        FrAve(i) = (Fr(i)+Fl(i+1))/2;
        FlAve(i) = Fl(i); 
    
        
    elseif i==N+1
        FrAve(i) = Fr(i);
        FlAve(i) = (Fr(i-1)+Fl(i))/2;
    else
  
      FrAve(i) = (Fr(i)+Fl(i+1))/2+jump(i+1);
      FlAve(i) = (Fr(i-1)+Fl(i))/2+jump(i);
      

    end

end

else
   assert(0) 
end
 
 if(strcmp(eqn,'solution')==1 || strcmp(eqn,'residual')==1)

     FI = (FrAve-FlAve)./h-obj.source;


elseif(strcmp(eqn,'error')==1)

  FI = (FrAve-FlAve)./h-obj.errorSource;


end

 
end






function [ FI ] = computeburgersviscfluxintegralb( obj,Z,eqn )
%COMPUTEFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here

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


Fr = zeros(N+2,1);
Fl = zeros(N+2,1);
FrAve = zeros(N+2,1);
FlAve = zeros(N+2,1);


if(obj.bcLeftType == 'D' && obj.bcRightType == 'D')
     
     ur = zeros(N+2,1);
        upr = zeros(N+2,1);
        ul = zeros(N+2,1);
        upl = zeros(N+2,1);
        utilder = zeros(N+2,1);
        utildel = zeros(N+2,1);
         utilderp = zeros(N+2,1);
        utildelp = zeros(N+2,1);
        
        nonlinearerror = 1;
         if(nonlinearerror && strcmp(eqn,'error')==1)
%              obj.computerespseudo();%
%              Zu = obj.unstructuredrecon(obj.convSoln,obj.qOrder,'error');
                      Zu = obj.convSolnRecon;%obj.unstructuredrecon(obj.convSoln,obj.qOrder,'error');
                      uorder = obj.qOrder;
%              obj.computeerrorpseudo();%
%                 obj.convSoln
%                 error('1')
 Zu;
% error('1')
         end
                
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
%    Fr(i) = -1*(ur(i)^2/2-upr(i)); %factor the flux function?
%    Fl(i) = -1*(ul(i)^2/2-upl(i));
     

    if(nonlinearerror && strcmp(eqn,'error')==1)
        for k = 1:uorder
        utilder(i) = utilder(i)+Zu(k,i)*(h(i)/2)^(k-1);
        utildel(i) = utildel(i)+Zu(k,i)*(-h(i)/2)^(k-1);
        end
        for k = 1:uorder-1
        utilderp(i) = utilderp(i) + k*Zu(k+1,i)*(h(i)/2)^(k-1);
        utildelp(i) = utildelp(i) + k*Zu(k+1,i)*(-h(i)/2)^(k-1);
        end


     Fr(i) = -1*((ur(i)+utilder(i))^2/2-(upr(i)+utilderp(i)) -( utilder(i)^2/2-utilderp(i) ) ); %factor the flux function?
     Fl(i) = -1*((ul(i)+utildel(i))^2/2-(upl(i)+utildelp(i)) -( utildel(i)^2/2-utildelp(i) ) );    


%   
% 
%       Fr(i) = Fr(i) - ur(i) * utilder(i);
%       Fl(i) = Fl(i) - ul(i) * utildel(i);
%       
    end

     if(p==2)
        jump(i) = (.2/((h(i)+h(i-1))/2))*(ul(i)-ur(i-1)) ;
      end

end

jump(2) = 0;
for i=2:N+1


    if i==2
        FrAve(i) = (Fr(i)+Fl(i+1))/2;
        FlAve(i) = Fl(i); 
    
        
    elseif i==N+1
        FrAve(i) = Fr(i);
        FlAve(i) = (Fr(i-1)+Fl(i))/2;
    else
  
      FrAve(i) = (Fr(i)+Fl(i+1))/2+jump(i+1);
      FlAve(i) = (Fr(i-1)+Fl(i))/2+jump(i);
      

    end

end

else
   assert(0) 
end
 
 if(strcmp(eqn,'solution')==1 || strcmp(eqn,'residual')==1)

     FI = (FrAve-FlAve)./h-obj.source;


elseif(strcmp(eqn,'error')==1)

  FI = (FrAve-FlAve)./h-obj.errorSource;


end

% obj.bcLeftVal
% obj.bcRightVal
% error('1')
 
end


