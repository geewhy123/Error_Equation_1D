function [ err ] = reconplot( obj,Z )
%RECONPLOT Summary of this function goes here
%   Detailed explanation goes here
% Z
% error('1')
err = 0;
x = obj.cellCentroids;
h = obj.cellWidths;
p = obj.pOrder;
N = obj.nCells;

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
ye = sin(pi*xx);
%  ye = (xx-0.5).^4;
err = max(err,max(abs(yy-ye)));
plot(xx,ye)
%     plot(xx,sin(pi*xx))
end


end

