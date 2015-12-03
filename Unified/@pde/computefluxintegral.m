function [ FI ] = computefluxintegral( obj,Z,eqn )
    if(strcmp(obj.physics,'Poisson')==1 || strcmp(obj.physics,'Advection')==1 || strcmp(obj.physics,'BurgersVisc')==1|| strcmp(obj.physics,'LinearSystem')==1||strcmp(obj.physics,'Burgers')==1 )    
        h = obj.cellWidths;
        N = obj.nCells;
        F = zeros(N+2,1);
        FI = zeros(N+2,1);
%        [err]= obj.reconplot(Z,'solution');
%         obj.curTime
% err
%         error('1')
        if(strcmp(eqn,'solution')==1 || strcmp(eqn,'residual')==1)
            f = obj.source;
%             if(strcmp(obj.goal,'TimeAccurate')==1)
%                 
%                 f = obj.computeTimeDepSource();
%             end
        elseif(strcmp(eqn,'error')==1)
            f = obj.errorSource;
        end
        
        for i = 1:N+1
            if(strcmp(obj.physics,'Poisson')==1)
                F(i) = computepoissonflux(obj,Z(:,i),Z(:,i+1),eqn,i);
%                 F
%                 error('1')
            elseif(strcmp(obj.physics,'Advection')==1)
                F(i) = computeadvectionflux(obj,Z(:,i),Z(:,i+1),eqn,i);
            elseif(strcmp(obj.physics,'BurgersVisc')==1)
                F(i) = computeburgersviscflux(obj,Z(:,i),Z(:,i+1),eqn,i);
            elseif(strcmp(obj.physics,'LinearSystem')==1)
              
               
%                Z
%                error('1')
               
               [ff] = computelinearsystemflux(obj,Z(:,i),Z(:,i+1),eqn,i);
               F1(i) = ff(1);
               F2(i) = ff(2);
            elseif(strcmp(obj.physics,'Burgers')==1)
                 F(i) = computeburgersflux(obj,Z(:,i),Z(:,i+1),eqn,i);
            end
            
            if(i==1)
                if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')
                    if(strcmp(obj.physics,'Poisson')==1)
                        F(i) = computepoissonflux(obj,Z(:,N+1),Z(:,i+1),eqn,i);
                    elseif(strcmp(obj.physics,'Advection')==1)
                        F(i) = computeadvectionflux(obj,Z(:,N+1),Z(:,i+1),eqn,i);
                    elseif(strcmp(obj.physics,'BurgersVisc')==1)
                        F(i) = computeburgersviscflux(obj,Z(:,N+1),Z(:,i+1),eqn,i);
                    elseif(strcmp(obj.physics,'LinearSystem')==1)
                        %                         F(i) = computelinearsystemflux(obj,Z(:,N+1),Z(:,i+1),eqn,i);
                        [ff] = computelinearsystemflux(obj,Z(:,N+1),Z(:,i+1),eqn,i);
                        F1(i) = ff(1);
                        F2(i) = ff(2);
                    elseif(strcmp(obj.physics,'Burgers')==1)
                        F(i) = computeburgersflux(obj,Z(:,N+1),Z(:,i+1),eqn,i);
                    end
                end
            elseif(i==N+1)
                if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')
                    F(i) = F(1);
                end
            end
        
            if(i>1)
                FI(i) = (F(i)-F(i-1))/h(i)-f(i);
                if(strcmp(obj.physics,'LinearSystem')==1)
                    FI1(i)=(F1(i)-F1(i-1))/h(i)-f(i);
                    FI2(i)=(F2(i)-F2(i-1))/h(i)-f(i);
                    FI(i,1) = FI1(i) ;
                    FI(i,2) = FI2(i) ;
                end
            end
        end
%     F
%     error('1')
        if(strcmp(eqn,'error')==1)
%              FI
%              if(isnan(FI(2)))
%                  Z
%                 error('1') 
%              end
%             Z
%             error('1')
%          error('1')
        end
      
    elseif(strcmp(obj.physics,'BurgersMod')==1)
        FI=computeburgersmodfluxintegral(obj,Z,eqn);
 
    else
        assert(0)
    end

% F1
% F2
end

function [ F ] = computepoissonflux( obj,left,right,eqn,i )
%COMPUTEFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here

    h = obj.cellWidths;
    N = obj.nCells;
    alpha = obj.jump(2);
    alphaL = obj.jump(1);
    alphaR = obj.jump(3);
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


% %     if(~isempty(obj.refinecells))
% %         if(i==1)
% %             p = obj.hOrder;
% %         elseif(i==N+1)
% %             p = obj.hOrder;
% %         else
% %             n = length(obj.refinecells);
% %             if(i < obj.refinecells(n/2)) 
% %                 pl = obj.hOrder;
% %                 pr = obj.hOrder;
% %             elseif(i==obj.refinecells(n/2))
% %                 pl = obj.hOrder;
% %                 pr = obj.pOrder;
% %             elseif(i== obj.refinecells(n/2+1)-1)
% %                 pl = obj.pOrder;
% %                 pr = obj.hOrder;
% %             elseif(i> obj.refinecells(n/2+1)-1)
% %                 pl = obj.hOrder;
% %                 pr = obj.hOrder;
% %             else
% %                 pl = obj.pOrder;
% %                 pr = obj.pOrder;
% %             end
% %             
% %         end
% %     else
    %     if(strcmp(eqn,'solution')==1)
    %     p = obj.pOrder;
        pr = p;%obj.pOrder;
        pl = p;%obj.pOrder;
    %     elseif(strcmp(eqn,'residual')==1)
    %     elseif(strcmp(eqn,'error')==1)
    %     end
% %     end       
 

    if(i==1 && obj.bcLeftType == 'D' && obj.bcRightType == 'D')
        F = 0;
        for k = 1:p-1
            F = F + k*right(k+1)*(-h(i+1)/2)^(k-1);
        end
          
        if(strcmp(obj.bchandle,'Jump')==1 )
            uu = 0;
            for k = 1:p
                uu = uu + right(k)*(-h(i+1)/2)^(k-1);
            end 
            F = F+alphaL*(uu-obj.bcLeftVal)/(h(i+1)/2);
            if(abs(uu-obj.bcLeftVal) > 1e-5)
%             [uu obj.bcLeftVal]
        %     error('1')
            end
        end
        %boundary jump
    
        return;
        
    elseif(i==N+1 && obj.bcLeftType == 'D' && obj.bcRightType == 'D')
        F = 0;
        for k = 1:p-1
            F = F + k*left(k+1)*(h(i)/2)^(k-1);
        end
        
        if(strcmp(obj.bchandle,'Jump')==1 )
            uu = 0;
            for k = 1:p
                uu = uu + left(k)*(h(i)/2)^(k-1);
            end 
            F = F-alphaR*(uu-obj.bcRightVal)/(h(i)/2);
            
%             error('1')
        end
        %boundary jump
        
        return;
        
    elseif(obj.bcLeftType=='P' && obj.bcRightType == 'P')
        Fl = 0;
        Fr = 0;
        
        for k = 1:p-1
            Fr = Fr + k*right(k+1)*(-h(i+1)/2)^(k-1);
        end
        for k = 1:p-1
            Fl = Fl + k*left(k+1)*(h(i)/2)^(k-1);
        end
        F = 0.5*(Fr+Fl);
    
        if(p==2 )
            ul1 = right(1)+right(2)*(-h(i+1)/2);
            ul2 = left(1) + left(2)*(h(i)/2);
            jump = (alpha/((h(i+1)+h(i))/2))*(ul1-ul2) ;
            F = F+jump;
                   
% %                       elseif(p==3)
% %           ul1 = right(1)+right(2)*(-h(i+1)/2)+right(3)*(-h(i+1)/2)^2;
% %             ul2 = left(1) + left(2)*(h(i)/2)+left(3)*(h(i)/2)^2;
% %                   jump = (alpha/((h(i+1)+h(i))/2))*(ul1-ul2) ;
% %                    F = F+jump;

        elseif(p==4)
            ul1 = right(1)+right(2)*(-h(i+1)/2)+right(3)*(-h(i+1)/2)^2+right(4)*(-h(i+1)/2)^3;
            ul2 = left(1) + left(2)*(h(i)/2)+left(3)*(h(i)/2)^2+left(4)*(h(i)/2)^3;
            jump = (1*alpha/((h(i+1)+h(i))/2))*(ul1-ul2) ;
            F = F+jump;
         
        end
        return;

    else
        Fl = 0;
        Fr = 0;
    
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
            jump = (alpha/((h(i+1)+h(i))/2))*(ul1-ul2) ;
%          end

            F = F+jump;
      
% %             elseif(pr==3 && pl==3)
% %           ul1 = right(1)+right(2)*(-h(i+1)/2)+right(3)*(-h(i+1)/2)^2;
% %             ul2 = left(1) + left(2)*(h(i)/2)+left(3)*(h(i)/2)^2;
% %                   jump = (alpha/((h(i+1)+h(i))/2))*(ul1-ul2) ;
% %                    F = F+jump;
         
      
        elseif(pr==4 && pl==4)
            ul1 = right(1)+right(2)*(-h(i+1)/2)+right(3)*(-h(i+1)/2)^2+right(4)*(-h(i+1)/2)^3;
            ul2 = left(1) + left(2)*(h(i)/2)+left(3)*(h(i)/2)^2+left(4)*(h(i)/2)^3;
            jump = (1*alpha/((h(i+1)+h(i))/2))*(ul1-ul2) ;
            F = F+jump;

    
        end
%     end
    
        return;
    
    end
 
end

function [ F ] = computeadvectionflux( obj,left,right,eqn,i  )
%COMPUTEFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here

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


    if(~isempty(obj.refinecells))
        if(i==1)
            p = obj.hOrder;
        elseif(i==N+1)
            p = obj.hOrder;
        else
            n = length(obj.refinecells);
            if(i < obj.refinecells(n/2)) 
                pl = obj.hOrder;
                pr = obj.hOrder;
            elseif(i==obj.refinecells(n/2))
                pl = obj.hOrder;
                pr = obj.pOrder;
            elseif(i== obj.refinecells(n/2+1)-1)
                pl = obj.pOrder;
                pr = obj.hOrder;
            elseif(i> obj.refinecells(n/2+1)-1)
                pl = obj.hOrder;
                pr = obj.hOrder;
            else
                pl = obj.pOrder;
                pr = obj.pOrder;
            end
            
        end
    else
    %     if(strcmp(eqn,'solution')==1)
    %     p = obj.pOrder;
        pr = p;%obj.pOrder;
        pl = p;%obj.pOrder;
    %     elseif(strcmp(eqn,'residual')==1)
    %     elseif(strcmp(eqn,'error')==1)
    %     end
    end       
 

    if(i==1 && obj.bcLeftType == 'F')% && obj.bcRightType == 'D')
        F = 0;
        for k = 1:p
            F = F + right(k)*(-h(i+1)/2)^(k-1);
        end
        return;
        
    elseif(i==N+1 && obj.bcRightType == 'D')% && obj.bcRightType == 'D')
        F = obj.bcRightVal;   
        return;
        
    elseif(obj.bcLeftType=='P' && obj.bcRightType == 'P')
        Fl = 0;
        Fr = 0;
    
        for k = 1:p
            Fr = Fr + right(k)*(-h(i+1)/2)^(k-1);
        end
        for k = 1:p
            Fl = Fl + left(k)*(h(i)/2)^(k-1);
        end
        F = Fr;%0.5*(Fr+Fl);
     
        if(p==2)
            ul1 = right(1)+right(2)*(-h(i+1)/2);
            ul2 = left(1) + left(2)*(h(i)/2);
            jump = (.0/((h(i+1)+h(i))/2))*(ul1-ul2) ;
            F = F+jump;
        end
%           if(i==5)
%             [left right]
%         [i Ul Ur]
%         error('1')
%         end
% [i Fl Fr]
% error('1')
        return;

    else
        Fl = 0;
        Fr = 0;
        
        for k = 1:pr
            Fr = Fr + right(k)*(-h(i+1)/2)^(k-1);
        end
        for k = 1:pl
            Fl = Fl + left(k)*(h(i)/2)^(k-1);
        end
        F =Fr;% 0.5*(Fr+Fl);
     
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
            jump = (.0/((h(i+1)+h(i))/2))*(ul1-ul2) ;
%          end

            F = F+jump;
        end
%     end
   
        return;
    
    end
    
 
end

function [ F ] = computeburgersflux( obj,left,right,eqn,i  )
%COMPUTEFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here

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
        
    if(obj.bcLeftType=='P' && obj.bcRightType == 'P')
        Ul = 0;
        Ur = 0;
        for k = 1:p
            Ur = Ur + right(k)*(-h(i+1)/2)^(k-1);
        end
        for k = 1:p
            Ul = Ul + left(k)*(h(i)/2)^(k-1);
        end
        Fr = -Ur^2/2;
        Fl = -Ul^2/2;
        dDot = (Ul+Ur)/2;
        F = 0.5*(Fr+Fl)-0.5*abs(dDot)*(Ul-Ur);

        return;

    else
        Fl = 0;
        Fr = 0;
        
        for k = 1:pr
            Fr = Fr + right(k)*(-h(i+1)/2)^(k-1);
        end
        for k = 1:pl
            Fl = Fl + left(k)*(h(i)/2)^(k-1);
        end
        F =Fl;%0.5*(Fr+Fl);%+abs(ddot)*(Ul-Ur);
     
        return;
    
    end
    
 
end


function [ F ] = computeburgersviscflux( obj,left,right,eqn,i  )
%COMPUTEFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here

    h = obj.cellWidths;
    N = obj.nCells;
    alpha = obj.jump(2);
    alphaL = obj.jump(1);
    alphaR = obj.jump(3);
    if(strcmp(eqn,'solution')==1)
        p = obj.pOrder;
    elseif(strcmp(eqn,'error')==1)
        p = obj.qOrder;
    elseif(strcmp(eqn,'residual')==1)
        p = obj.rOrder;
    else
        assert(0); 
    end


    if(~isempty(obj.refinecells))
        if(i==1)
            p = obj.hOrder;
        elseif(i==N+1)
            p = obj.hOrder;
        else
            n = length(obj.refinecells);
            if(i < obj.refinecells(n/2)) 
                pl = obj.hOrder;
                pr = obj.hOrder;
            elseif(i==obj.refinecells(n/2))
                pl = obj.hOrder;
                pr = obj.pOrder;
            elseif(i== obj.refinecells(n/2+1)-1)
                pl = obj.pOrder;
                pr = obj.hOrder;
            elseif(i> obj.refinecells(n/2+1)-1)
                pl = obj.hOrder;
                pr = obj.hOrder;
            else
                pl = obj.pOrder;
                pr = obj.pOrder;
            end
            
        end
    else
%     if(strcmp(eqn,'solution')==1)
%     p = obj.pOrder;
        pr = p;%obj.pOrder;
        pl = p;%obj.pOrder;
%     elseif(strcmp(eqn,'residual')==1)
%     elseif(strcmp(eqn,'error')==1)
%     end
    end       
 

    nonlinearerror = ~strcmp(obj.NLError,'PrimalError');
    utilder = 0;
    utildel = 0;
    if(i==1 && obj.bcLeftType == 'D' && obj.bcRightType == 'D')
        F = 0;
        for k = 1:p-1
            F = F + k*right(k+1)*(-h(i+1)/2)^(k-1);
        end
        U=0;
        for k = 1:p
            U = U+(right(k)*(-h(i+1)/2)^(k-1));
        end
        
        %
        if(strcmp(obj.NLError,'NewtonError')==1 && strcmp(eqn,'error')==1)
           F = F+U^2/2; 
        end
        %
        
        F = F-U^2/2;
   
        
     
        if(nonlinearerror && strcmp(eqn,'error')==1)
            Zu = obj.convSolnRecon;
            uorder = obj.qOrder;
%             error('1')
        end
        
        if(nonlinearerror && strcmp(eqn,'error')==1)
            for k = 1:uorder
    %         utildel = utildel+Zu(k,i)*(h(i)/2)^(k-1);
                utilder = utilder+Zu(k,i+1)*(-h(i+1)/2)^(k-1);
            end
            F = F - U * utilder;
    %         Fl = Fl - ul * utildel;
           
    
    %
%         if(strcmp(eqn,'error')==1)
%         fprintf('warning - linearized flux')
%         F = F + U^2/2;
%         end
        %
        end
        return;
        
    elseif(i==N+1 && obj.bcLeftType == 'D' && obj.bcRightType == 'D')
        F = 0;
        for k = 1:p-1
            F = F + k*left(k+1)*(h(i)/2)^(k-1);
        end
        U=0;
        for k = 1:p
            U = U+(left(k)*(h(i)/2)^(k-1)) ;
        end
        
        
        %
         
        if(strcmp(obj.NLError,'NewtonError')==1 && strcmp(eqn,'error')==1)
           F = F+U^2/2; 
        end
        
        %
        
        F = F-U^2/2;
   
        if(nonlinearerror && strcmp(eqn,'error')==1)
            Zu = obj.convSolnRecon;
            uorder = obj.qOrder;
        end
        
        if(nonlinearerror && strcmp(eqn,'error')==1)
            for k = 1:uorder
                utildel = utildel+Zu(k,i)*(h(i)/2)^(k-1);
%         utilder = utilder+Zu(k,i+1)*(-h(i+1)/2)^(k-1);
            end
%         F = F - ur * utilder;
            F = F - U * utildel;
   %
%         if(strcmp(eqn,'error')==1)
%         fprintf('warning - linearized flux')
%         F = F + U^2/2;
%         end
        %
        end
        return;
        
    elseif(obj.bcLeftType=='P' && obj.bcRightType == 'P')
        Fl = 0;
        Fr = 0;
        
        ur=0;
        for k = 1:p-1
            Fr = Fr + k*right(k+1)*(-h(i+1)/2)^(k-1);
        end
        for k = 1:p
            ur = ur+(right(k)*(-h(i+1)/2)^(k-1));
        end
        Fr = Fr-ur^2/2;
    
        ul=0;
        for k = 1:p-1
            Fl = Fl + k*left(k+1)*(h(i)/2)^(k-1);
        end
        for k = 1:p
            ul= ul+ left(k)*(h(i)/2)^(k-1);
        end
        Fl = Fl-ul^2/2;
     
        if(nonlinearerror && strcmp(eqn,'error')==1)
            Zu = obj.convSolnRecon;
            uorder = obj.qOrder;
        end
        if(nonlinearerror && strcmp(eqn,'error')==1)
            for k = 1:uorder
                utildel = utildel+Zu(k,i)*(h(i)/2)^(k-1);
                utilder = utilder+Zu(k,i+1)*(-h(i+1)/2)^(k-1);
            end
            Fr = Fr - ur * utilder;
            Fl = Fl - ul * utildel;

        end
  
        
        %
%         if(strcmp(eqn,'error')==1)
%         fprintf('warning - linearized flux')
%         Fr = Fr + ur^2/2;
%         Fl = Fl + ul^2/2;

%         end
        %
        F = 0.5*(Fr+Fl);


        if(p==2)
            ul1 = right(1)+right(2)*(-h(i+1)/2);
            ul2 = left(1) + left(2)*(h(i)/2);
            jump = (alpha/((h(i+1)+h(i))/2))*(ul1-ul2) ;
            F = F+jump;
        end
        return;
        
    else
        Fl = 0;
        Fr = 0;
        for k = 1:pr-1
            Fr = Fr + k*right(k+1)*(-h(i+1)/2)^(k-1);
        end
        ur=0;
        for k = 1:pr
            ur=ur+(right(k)*(-h(i+1)/2)^(k-1)); 
        end
        
                %
        if(strcmp(obj.NLError,'NewtonError')==1 && strcmp(eqn,'error')==1)
           Fr = Fr+ur^2/2; 
        end
        %
        
        Fr = Fr-ur^2/2;
        for k = 1:pl-1
            Fl = Fl + k*left(k+1)*(h(i)/2)^(k-1);
        end
        ul=0;
        for k = 1:pl
            ul = ul+(left(k)*(h(i)/2)^(k-1));
        end
        
        
                %
        if(strcmp(obj.NLError,'NewtonError')==1 && strcmp(eqn,'error')==1)
           Fl = Fl+ul^2/2; 
        end
        %
        
        Fl = Fl-ul^2/2;
      
      

        if(nonlinearerror && strcmp(eqn,'error')==1)
            Zu = obj.convSolnRecon;
            uorder = obj.qOrder;
        end
        if(nonlinearerror && strcmp(eqn,'error')==1)
            for k = 1:uorder
                utildel = utildel+Zu(k,i)*(h(i)/2)^(k-1);
                utilder = utilder+Zu(k,i+1)*(-h(i+1)/2)^(k-1);
            end
            Fr = Fr - ur * utilder;
            Fl = Fl - ul * utildel;

        end
       %
%         if(strcmp(eqn,'error')==1)
%         fprintf('warning - linearized flux')
%         Fr = Fr + ur^2/2;
%         Fl = Fl + ul^2/2;
% 
%         end
        %
        F = 0.5*(Fr+Fl);

        if(pr==2 && pl == 2)
%          error('1')
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
            jump = (alpha/((h(i+1)+h(i))/2))*(ul1-ul2) ;
%          end

            F = F+jump;
        end
%     end
   
        return;
    
    end
   
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




% function [ FI ] = computeburgersviscfluxintegralb( obj,Z,eqn )
% %COMPUTEFLUXINTEGRAL Summary of this function goes here
% %   Detailed explanation goes here
% 
% if(strcmp(eqn,'solution')==1)
%     p = obj.pOrder;
% elseif(strcmp(eqn,'error')==1)
%     p = obj.qOrder;
% elseif(strcmp(eqn,'residual')==1)
%     p = obj.rOrder;
% else
%    assert(0); 
% end
% 
% x = obj.cellCentroids;
% h = obj.cellWidths;
% N = obj.nCells;
% 
% 
% Fr = zeros(N+2,1);
% Fl = zeros(N+2,1);
% FrAve = zeros(N+2,1);
% FlAve = zeros(N+2,1);
% 
% 
% if(obj.bcLeftType == 'D' && obj.bcRightType == 'D')
%      
%      ur = zeros(N+2,1);
%         upr = zeros(N+2,1);
%         ul = zeros(N+2,1);
%         upl = zeros(N+2,1);
%         utilder = zeros(N+2,1);
%         utildel = zeros(N+2,1);
%          utilderp = zeros(N+2,1);
%         utildelp = zeros(N+2,1);
%         
%         nonlinearerror = 1;
%          if(nonlinearerror && strcmp(eqn,'error')==1)
% %              obj.computerespseudo();%
% %              Zu = obj.unstructuredrecon(obj.convSoln,obj.qOrder,'error');
%                       Zu = obj.convSolnRecon;%obj.unstructuredrecon(obj.convSoln,obj.qOrder,'error');
%                       uorder = obj.qOrder;
% %              obj.computeerrorpseudo();%
% %                 obj.convSoln
% %                 error('1')
%  Zu;
% % error('1')
%          end
%                 
% jump = zeros(N+2,1);
% for i=2:N+1
% 
%     for k = 1:p
%         ur(i) = ur(i)+Z(k,i)*(h(i)/2)^(k-1);
%         ul(i) = ul(i)+Z(k,i)*(-h(i)/2)^(k-1);
%     end
%     for k = 1:p-1
%         upr(i) = upr(i)+k*Z(k+1,i)*(h(i)/2)^(k-1);
%         upl(i) = upl(i)+k*Z(k+1,i)*(-h(i)/2)^(k-1);
%     end
% 
% %         assert(ur(i)-1 < 0);
% %    Fr(i) = -1*(ur(i)^2/2-upr(i)); %factor the flux function?
% %    Fl(i) = -1*(ul(i)^2/2-upl(i));
%      
% 
%     if(nonlinearerror && strcmp(eqn,'error')==1)
%         for k = 1:uorder
%         utilder(i) = utilder(i)+Zu(k,i)*(h(i)/2)^(k-1);
%         utildel(i) = utildel(i)+Zu(k,i)*(-h(i)/2)^(k-1);
%         end
%         for k = 1:uorder-1
%         utilderp(i) = utilderp(i) + k*Zu(k+1,i)*(h(i)/2)^(k-1);
%         utildelp(i) = utildelp(i) + k*Zu(k+1,i)*(-h(i)/2)^(k-1);
%         end
% 
% 
%      Fr(i) = -1*((ur(i)+utilder(i))^2/2-(upr(i)+utilderp(i)) -( utilder(i)^2/2-utilderp(i) ) ); %factor the flux function?
%      Fl(i) = -1*((ul(i)+utildel(i))^2/2-(upl(i)+utildelp(i)) -( utildel(i)^2/2-utildelp(i) ) );    
% 
% 
% %   
% % 
% %       Fr(i) = Fr(i) - ur(i) * utilder(i);
% %       Fl(i) = Fl(i) - ul(i) * utildel(i);
% %       
%     end
% 
%      if(p==2)
%         jump(i) = (.2/((h(i)+h(i-1))/2))*(ul(i)-ur(i-1)) ;
%       end
% 
% end
% 
% jump(2) = 0;
% for i=2:N+1
% 
% 
%     if i==2
%         FrAve(i) = (Fr(i)+Fl(i+1))/2;
%         FlAve(i) = Fl(i); 
%     
%         
%     elseif i==N+1
%         FrAve(i) = Fr(i);
%         FlAve(i) = (Fr(i-1)+Fl(i))/2;
%     else
%   
%       FrAve(i) = (Fr(i)+Fl(i+1))/2+jump(i+1);
%       FlAve(i) = (Fr(i-1)+Fl(i))/2+jump(i);
%       
% 
%     end
% 
% end
% 
% else
%    assert(0) 
% end
%  
%  if(strcmp(eqn,'solution')==1 || strcmp(eqn,'residual')==1)
% 
%      FI = (FrAve-FlAve)./h-obj.source;
% 
% 
% elseif(strcmp(eqn,'error')==1)
% 
%   FI = (FrAve-FlAve)./h-obj.errorSource;
% 
% 
% end
% 
% % obj.bcLeftVal
% % obj.bcRightVal
% % error('1')
%  
% end

function [ F ] = computelinearsystemflux( obj,left,right,eqn,i  )
%COMPUTEFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here

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


    if(~isempty(obj.refinecells))
        if(i==1)
            p = obj.hOrder;
        elseif(i==N+1)
            p = obj.hOrder;
        else
            n = length(obj.refinecells);
            if(i < obj.refinecells(n/2)) 
                pl = obj.hOrder;
                pr = obj.hOrder;
            elseif(i==obj.refinecells(n/2))
                pl = obj.hOrder;
                pr = obj.pOrder;
            elseif(i== obj.refinecells(n/2+1)-1)
                pl = obj.pOrder;
                pr = obj.hOrder;
            elseif(i> obj.refinecells(n/2+1)-1)
                pl = obj.hOrder;
                pr = obj.hOrder;
            else
                pl = obj.pOrder;
                pr = obj.pOrder;
            end
            
        end
    else
    %     if(strcmp(eqn,'solution')==1)
    %     p = obj.pOrder;
        pr = p;%obj.pOrder;
        pl = p;%obj.pOrder;
    %     elseif(strcmp(eqn,'residual')==1)
    %     elseif(strcmp(eqn,'error')==1)
    %     end
    end       
 

    if(i==1 && obj.bcLeftType == 'F')% && obj.bcRightType == 'D')
        F(1) = 0;
        F(2) = 0;
        for k = 1:p
            F(1) = F(1) + right(k)*(-h(i+1)/2)^(k-1);
            F(2) = F(2) + right(k+p)*(-h(i+1)/2)^(k-1);
        end
        return;
        
    elseif(i==N+1 && obj.bcRightType == 'D')% && obj.bcRightType == 'D')
        F(1) = obj.bcRightVal(1);
        F(2) = obj.bcRightVal(2);
        return;
        
    elseif(obj.bcLeftType=='P' && obj.bcRightType == 'P')
        Fl = 0;
        Fr = 0;
    
        for k = 1:p
            Fr = Fr + right(k)*(-h(i+1)/2)^(k-1);
        end
        for k = 1:p
            Fl = Fl + left(k)*(h(i)/2)^(k-1);
        end
        F = Fr;%0.5*(Fr+Fl);
     
        if(p==2)
            ul1 = right(1)+right(2)*(-h(i+1)/2);
            ul2 = left(1) + left(2)*(h(i)/2);
            jump = (.0/((h(i+1)+h(i))/2))*(ul1-ul2) ;
            F = F+jump;
        end
        return;

    else
       
        Fl(1) = 0;
        Fl(2) = 0;
        Fr(1) = 0;
        Fr(2) = 0;
        
        for k = 1:pr
            Fr(2) = Fr(1) + right(k)*(-h(i+1)/2)^(k-1);
            Fr(1) = Fr(2) + right(k+p)*(-h(i+1)/2)^(k-1);
        end
        for k = 1:pl
            Fl(2) = Fl(1) + left(k)*(h(i)/2)^(k-1);
            Fl(1) = Fl(2) + left(k+p)*(h(i)/2)^(k-1);
        end
        F =Fr;% 0.5*(Fr+Fl);
     
%         if(pr==2 && pl == 2)
% %          if(i==1)
% %              if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')
% %                 ul1 = right(1)+right(2)*(-h(i+1)/2);
% %                 ul2 = left(1)+ left(2)*(h(N+1)/2);
% %                 jump = (.2/((h(i)+h(N+1))/2))*(ul1-ul2) ;
% %              end
% %          elseif(i==N+1)
% %                 ul1 = right(1)+right(2)*(-h(i)/2);
% %                 ul2 = left(1)+ left(2)*(h()/2);
% %                 jump = (.2/((h(i)+h(N+1))/2))*(ul1-ul2) ;
% %              if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')
% %                  
% %              end
% %          
% %          else
%             ul1 = right(1)+right(2)*(-h(i+1)/2);
%             ul2 = left(1) + left(2)*(h(i)/2);
%             jump = (.0/((h(i+1)+h(i))/2))*(ul1-ul2) ;
% %          end
% 
%             F = F+jump;
%         end
%     end
    
        return;
    
    end
    
 
end


