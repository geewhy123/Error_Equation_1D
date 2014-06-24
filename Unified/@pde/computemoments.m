function  computemoments( obj)
%COMPUTEMOMENTS Summary of this function goes here
%   Detailed explanation goes here


 x = obj.cellCentroids;
 h = obj.cellWidths;
 N = obj.nCells;
 
 
m = 6;%max(max(obj.pOrder,max(obj.qOrder,obj.rOrder)),obj.hOrder);
 A = NaN*ones(N+2,m);
 for i = 2:N+1
    for j = 1 : m+1
       
A(i,j) = (1/h(i))*( ((h(i)/2))^(j)/(j) -((-h(i)/2))^(j)/(j) ); 
        
    end
     
     
 end


 obj.moments = A;


end

