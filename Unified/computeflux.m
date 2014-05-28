function [left,right ] = computeflux(Z,h,N,p,phys,equation,obj )
%COMPUTEFLUX Summary of this function goes here
%   Detailed explanation goes here
% error('1')


% global dir
if (obj.bcLeftType == 'D' && obj.bcRightType == 'D')

        
        
        
left  = NaN*ones(N+2,1);
right = NaN*ones(N+2,1);
up = zeros(N+1,1);

for i = 2:N+2
 if (i==2)
        y = Z(:,i);
%         yl = -y;%zeros(size(Z(:,i)));
       % yl = Z(:,N+1);
       
 elseif(i==N+2)

%         y = Z(:,2);

        yl = Z(:,i-1);
%         y = -yl;%zeros(size(Z(:,i)));
 else
        y = Z(:,i);
        yl = Z(:,i-1);
    
end    

switch p
    case 2
        
        if(i==2)
        yl = 0;    
        upl1 = y(2);
        upl2 = upl1;
        ul1 = 0;
        ul2 = 0;
        elseif(i==N+2)
            y=0;
        upl2 = yl(2);
        upl1 = upl2;
        ul1 = 0;
        ul2 = 0;
        else
        
        upl1 = y(2);
        upl2 = yl(2);
        
        
%         upl1
        
        ul1 = y(1)+y(2)*(-h(i)/2);
        ul2 = yl(1) + yl(2)*(h(i-1)/2);
        end
        upl = (upl1+upl2)/2 + (.2/((h(i)+h(i-1))/2))*(ul1-ul2);

% if(i==4)
%     upl
%     error('1')
% end
        

    case 3
       
  if(i==2)
        yl = 0;    
        upl1 =  y(2)+2*y(3)*-h(i)/2;
        upl2 = upl1;
        ul1 = 0;
        ul2 = 0;
        elseif(i==N+2)
            y=0;
        upl2 = yl(2)+2*yl(3)*h(i-1)/2;
        upl1 = upl2;
        ul1 = 0;
        ul2 = 0;
        else

      
        upl1 = y(2)+2*y(3)*(-h(i)/2);
        upl2 = yl(2)+2*yl(3)*h(i-1)/2;
        ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 ;
        ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 ;
  end
        upl = (upl1+upl2)/2;% +(.2/((h(i)+h(i+1))/2))*(ul1-ul2);
        

    case 4
         if(i==2)
        yl = 0;    
        
        upl1 =  y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
        upl2 = upl1;
        ul1 = 0;
        ul2 = 0;
        elseif(i==N+2)
            y=0;
        upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2;
        upl1 = upl2;
        ul1 = 0;
        ul2 = 0;
        else

        upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
        upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2;

        ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3;
        ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 + yl(4)*(h(i-1)/2)^3;
         end
        upl = (upl1+upl2)/2;% +(.2/((h(i)+h(i+1))/2))*(ul1-ul2);
             

    case 5
       if(i==2)
        yl = 0;    
        upl1 =  y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2+4*y(5)*(-h(i)/2)^3;
        upl2 = upl1;
        ul1 = 0;
        ul2 = 0;
        elseif(i==N+2)
            y=0;
        upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2+4*yl(5)*(h(i-1)/2)^3;
        upl1 = upl2;
        ul1 = 0;
        ul2 = 0;
        else

        upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2   + 4*y(5)*(-h(i)/2)^3;
        upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2   + 4*yl(5)*(h(i-1)/2)^3;
        ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3+y(5)*(-h(i)/2)^4;
        ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 + yl(4)*(h(i-1)/2)^3 +yl(5)*(h(i-1)/2)^4; 
       end
        upl = (upl1+upl2)/2;


    case 6
        if(i==2)
        yl = 0;    
        upl1 =  y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2+4*y(5)*(-h(i)/2)^3+5*y(6)*(-h(i)/2)^4;
        upl2 = upl1;
        ul1 = 0;
        ul2 = 0;
        elseif(i==N+2)
            y=0;
        upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2+4*yl(5)*(h(i-1)/2)^3+5*yl(6)*(h(i-1)/2)^4;
        upl1 = upl2;
        ul1 = 0;
        ul2 = 0;
        else


        upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2   + 4*y(5)*(-h(i)/2)^3+ 5*y(6)*(-h(i)/2)^4;
        upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2   + 4*yl(5)*(h(i-1)/2)^3 + 5*yl(6)*(h(i-1)/2)^4; 
        ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3+y(5)*(-h(i)/2)^4 +y(6)*(-h(i)/2)^5;
        ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 + yl(4)*(h(i-1)/2)^3 +yl(5)*(h(i-1)/2)^4 + yl(6)*(h(i-1)/2)^5; 
        end
        upl = (upl1+upl2)/2;
           
        
    otherwise
        error('1')
end

switch phys
    case 'Poisson'
        up(i-1) = upl;
    case 'Advection'
       up(i-1) = ul1;

    case 'Burgers'



up(i-1) = ul2;
if((ul1+ul2)<0)
    up(i-1) = ul1;
end

end



end

left(2:N+1) = up(1:N);
right(2:N+1) = up(2:N+1);



    
else
    
    
left  = NaN*ones(N+2,1);
right = NaN*ones(N+2,1);
up = zeros(N+1,1);

for i = 2:N+2
 if (i==2)
        y = Z(:,i);
        yl = Z(:,N+1);
 elseif(i==N+2)
        y = Z(:,2);
        yl = Z(:,i-1);
 else
        y = Z(:,i);
        yl = Z(:,i-1);
    
end    

switch p
    case 2
          
        upl1 = y(2);
        upl2 = yl(2);
        ul1 = y(1)+y(2)*(-h(i)/2);
        ul2 = yl(1) + yl(2)*(h(i-1)/2);
        upl = (upl1+upl2)/2 + (.2/((h(i)+h(i-1))/2))*(ul1-ul2);


    case 3
       


      
        upl1 = y(2)+2*y(3)*(-h(i)/2);
        upl2 = yl(2)+2*yl(3)*h(i-1)/2;
        ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 ;
        ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 ;
        upl = (upl1+upl2)/2;% +(.2/((h(i)+h(i+1))/2))*(ul1-ul2);
        

    case 4
       

        upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
        upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2;

        ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3;
        ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 + yl(4)*(h(i-1)/2)^3;
        upl = (upl1+upl2)/2;% +(.2/((h(i)+h(i+1))/2))*(ul1-ul2);
             

    case 5
       

        upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2   + 4*y(5)*(-h(i)/2)^3;
        upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2   + 4*yl(5)*(h(i-1)/2)^3;
        ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3+y(5)*(-h(i)/2)^4;
        ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 + yl(4)*(h(i-1)/2)^3 +yl(5)*(h(i-1)/2)^4; 
        upl = (upl1+upl2)/2;


    case 6
      

        upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2   + 4*y(5)*(-h(i)/2)^3+ 5*y(6)*(-h(i)/2)^4;
        upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2   + 4*yl(5)*(h(i-1)/2)^3 + 5*yl(6)*(h(i-1)/2)^4; 
        ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3+y(5)*(-h(i)/2)^4 +y(6)*(-h(i)/2)^5;
        ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 + yl(4)*(h(i-1)/2)^3 +yl(5)*(h(i-1)/2)^4 + yl(6)*(h(i-1)/2)^5; 
        upl = (upl1+upl2)/2;
           
        
    otherwise
        error('1')
end

switch phys
    case 'Poisson'
        up(i-1) = upl;
    case 'Advection'
       up(i-1) = ul1;

    case 'Burgers'
%         tmp1 = ul1;
%         tmp2 = ul2;
%         ul1 = tmp2;
%         ul2 = tmp1;
        
        
%          up(i-1) = ul2;
%          if(ul1>ul2)
%             S = 0.5*(ul1+ul2);
%             if(S>0)
%                 up(i-1) = ul1;
%             else
%                up(i-1)=ul2; 
%             end
%                 
%          else
%              if(ul1>0)
%                  up(i-1) = ul1;
%              elseif(ul2<0)
%                  up(i-1) = ul2;
%              else
%                  up(i-1)=0;
%              end
%              
%          end



up(i-1) = ul2;
if((ul1+ul2)<0)
    up(i-1) = ul1;
end

% % if(strcmp(equation,'solution')==1 || strcmp(equation,'residual')==1)
% %             if(abs(ul1+ul2) < 1e-12 )
% %                up(i-1) = 0;
% %                 dir(i) = 0; 
% %             elseif((ul1+ul2)>0 )
% %                 up(i-1) = ul2;
% %                 dir(i) = 1;
% % 
% %  
% %             elseif( (ul1+ul2)<0 )
% %                 up(i-1) = ul1;
% %                 dir(i) = -1;
% %             else
% %                 error('5');
% %             end
% %     
% % elseif(strcmp(equation,'error')==1)
% %             if(dir(i)==0 )
% %                up(i-1) = 0;
% %                 dir(i) = 0; 
% %             elseif(dir(i)==1 )
% %                 up(i-1) = ul2;
% %                 dir(i) = 1;
% % 
% %  
% %             elseif( dir(i)==-1 )
% %                 up(i-1) = ul1;
% %                 dir(i) = -1;
% %             else
% %                 error('5');
% %             end
% %     
% %     
% % else   error('6')
% %     
% % end
            
%             if(abs(ul1+ul2) < 10^-6)
%                up(i-1) = 0; 
%             end
end



end

left(2:N+1) = up(1:N);
right(2:N+1) = up(2:N+1);



% for i = 2:N+1
%  if (i==2)
%         y = Z(:,i);
%         yr = Z(:,i+1);
%         yl = Z(:,N+1);
%  elseif(i== N+1)
%         y = Z(:,i);
%         yr = Z(:,2);
%         yl = Z(:,i-1);
%  else
%         y = Z(:,i);
%         yr = Z(:,i+1);
%         yl = Z(:,i-1);
%     
% end    
% 
% switch p
%     case 2
%         upr1 = y(2);
%         upr2 = yr(2);
%         ur2 = yr(1)+yr(2)*(-h(i+1)/2);
%         ur1 = y(1) + y(2)*(h(i)/2);
%         upr = (upr1+upr2)/2 + (.2/((h(i)+h(i+1))/2))*(ur2-ur1);
%         
%           
%         upl1 = y(2);
%         upl2 = yl(2);
%         ul1 = y(1)+y(2)*(-h(i)/2);
%         ul2 = yl(1) + yl(2)*(h(i-1)/2);
%         upl = (upl1+upl2)/2 + (.2/((h(i)+h(i-1))/2))*(ul1-ul2);
% 
% 
%     case 3
%         upr1 = y(2)+2*y(3)*h(i)/2; % u_i+1/2 using recon in i 
%         upr2 = yr(2)+2*yr(3)*(-h(i+1)/2);
%         ur2 = yr(1)+yr(2)*(-h(i+1)/2)+yr(3)*(-h(i+1)/2)^2;
%         ur1 = y(1) + y(2)*(h(i)/2)+y(3)*(h(i)/2)^2;
% 
% 
%         upr = (upr1+upr2)/2 ;%+(.2/((h(i)+h(i+1))/2))*(ur2-ur1);
%         upl1 = y(2)+2*y(3)*(-h(i)/2);
%         upl2 = yl(2)+2*yl(3)*h(i-1)/2;
%         ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 ;
%         ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 ;
%         upl = (upl1+upl2)/2;% +(.2/((h(i)+h(i+1))/2))*(ul1-ul2);
%         
% 
%     case 4
%         upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
%         upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2;
%         ur2 = yr(1)+yr(2)*(-h(i+1)/2)+yr(3)*(-h(i+1)/2)^2+yr(4)*(-h(i+1)/2)^3;
%         ur1 = y(1) + y(2)*(h(i)/2)+y(3)*(h(i)/2)^2 + y(4)*(h(i)/2)^3;
%         upr = (upr1+upr2)/2 ;%+(.2/((h(i)+h(i+1))/2))*(ur2-ur1);
% 
%         upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
%         upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2;
% 
%         ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3;
%         ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 + yl(4)*(h(i-1)/2)^3;
%         upl = (upl1+upl2)/2;% +(.2/((h(i)+h(i+1))/2))*(ul1-ul2);
%              
% 
%     case 5
%         upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2  + 4*y(5)*(h(i)/2)^3;
%         upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2   + 4*yr(5)*(-h(i+1)/2)^3;
%         ur2 = yr(1)+yr(2)*(-h(i+1)/2)+yr(3)*(-h(i+1)/2)^2+yr(4)*(-h(i+1)/2)^3+yr(5)*(-h(i+1)/2)^4;
%         ur1 = y(1) + y(2)*(h(i)/2)+y(3)*(h(i)/2)^2 + y(4)*(h(i)/2)^3 + y(5)*(h(i)/2)^4;
%         upr = (upr1+upr2)/2;
% 
%         upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2   + 4*y(5)*(-h(i)/2)^3;
%         upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2   + 4*yl(5)*(h(i-1)/2)^3;
%         ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3+y(5)*(-h(i)/2)^4;
%         ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 + yl(4)*(h(i-1)/2)^3 +yl(5)*(h(i-1)/2)^4; 
%         upl = (upl1+upl2)/2;
% 
% 
%     case 6
%         upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2  + 4*y(5)*(h(i)/2)^3+5*y(6)*(h(i)/2)^4;
%         upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2   + 4*yr(5)*(-h(i+1)/2)^3 + 5*yr(6)*(-h(i+1)/2)^4;
%         ur2 = yr(1)+yr(2)*(-h(i+1)/2)+yr(3)*(-h(i+1)/2)^2+yr(4)*(-h(i+1)/2)^3+yr(5)*(-h(i+1)/2)^4 +yr(6)*(-h(i+1)/2)^5;
%         ur1 = y(1) + y(2)*(h(i)/2)+y(3)*(h(i)/2)^2 + y(4)*(h(i)/2)^3 + y(5)*(h(i)/2)^4 +y(6)*(h(i)/2)^5;
%         upr = (upr1+upr2)/2;
% 
%         upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2   + 4*y(5)*(-h(i)/2)^3+ 5*y(6)*(-h(i)/2)^4;
%         upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2   + 4*yl(5)*(h(i-1)/2)^3 + 5*yl(6)*(h(i-1)/2)^4; 
%         ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3+y(5)*(-h(i)/2)^4 +y(6)*(-h(i)/2)^5;
%         ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 + yl(4)*(h(i-1)/2)^3 +yl(5)*(h(i-1)/2)^4 + yl(6)*(h(i-1)/2)^5; 
%         upl = (upl1+upl2)/2;
%            
%         
%     otherwise
%         error('1')
% end
% 
% switch phys
%     case 'Poisson'
%         left(i) = upl;
%         right(i) = upr;
%     case 'Advection'
%         left(i) = ul1;
%         right(i) = ur2;
%     case 'Burgers'
%         left(i) = ul2;
%         right(i) = ur1;
%     otherwise
%         error('5')
%     
% end
% 
% end
% up

% left
% right
% global xx
% size(xx)
% size(up)
% max(abs(up-sin(2*pi*xx)/(2*pi)))
% plot(xx,up-sin(2*pi*xx)/(2*pi),'o-')
% error('1')

end

end

