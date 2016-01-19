function [ err ] = reconplot( obj,Z ,eqn)
%RECONPLOT Summary of this function goes here
%   Detailed explanation goes here
% Z
% error('1')
err = 0;
errp=0;
x = obj.cellCentroids;
h = obj.cellWidths;

if(strcmp(eqn,'solution')==1)
p = max(obj.pOrder,obj.hOrder);
elseif(strcmp(eqn,'residual')==1)
p = obj.rOrder;
elseif(strcmp(eqn,'error')==1)
p = obj.qOrder;
elseif(strcmp(eqn,'average')==1)
p = 1;
end


% N = obj.nCells;
[~,N] = size(Z(:,2:end-1));
N

for i = 2:N+1
   xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,100);
   yy = 0;
   order = p;
   j=1;
   while(order >0 )
       yy = yy + Z(j,i)*(xx-x(i)).^(j-1);
    
       order = order-1 ;
      j=j+1;
      

   end
   
      
plot(xx,yy)
hold on
 ye=2*pi*cos(2*pi*xx)./(sin(2*pi*xx)+exp(1)^3);
 ye = sin(pi*xx)+1;
 ye = 1;
% ye = sin(pi*xx);
% ye = xx.^2.*(1-xx);
ye = 1-tanh(xx/2);

ye = sin(pi*xx);

ye = (1/(2*pi*2))*sin(pi*xx);
%  ye = (xx-0.5).^4;
err = max(err,max(abs(yy-ye)));
% plot(xx,ye)
%     plot(xx,sin(pi*xx))

end



% %derivative
% figure
% Z
% for i = 2:N+1
%    xx = linspace(x(i)-h(i)/2,x(i)+h(i)/2,100);
%    yyp = 0;
%    order = p;
%    j=1;
%    while(order >1 )
%        yyp = yyp + j*Z(j+1,i)*(xx-x(i)).^(j-1);
%     
%        order = order-1 ;
%       j=j+1;
%       
% 
%    end
%    
%       
% plot(xx,yyp)
% hold on
%  yep=2*pi*cos(2*pi*xx)./(sin(2*pi*xx)+exp(1)^3);
%  yep = sin(pi*xx)+1;
%  yep = 2*xx-3*xx.^2;
% % yep = sin(pi*xx);
% % yep = xx.^2.*(1-xx);
% % yep = 1-tanh(xx/2);
% %  yep = (xx-0.5).^4;
% errp = max(errp,max(abs(yyp-yep)));
%  plot(xx,yep)
% end
% 
% errp
% error('1')



end

