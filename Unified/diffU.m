    function [ dUdt] = diffU( U ,k)
%DIFFU Summary of this function goes here
%   Detailed explanation goes here

[m,n] = size(U);



% 6th order
%A = [1 1 1 1 1 1 1 1 1 ; -4 -3 -2 -1 0 1 2 3 4; 
%b = [0 1 0]';

%centered
for j = 1:9
    for i = 1:9
     A(i,j) = (j-5)^(i-1)/factorial(i-1);            
    end
end
                
    b = zeros(9,1);
    b(2) = 1;
cc = A\b;

for j = 1:9
    for i = 1:9
     Al(i,j) = (j-1)^(i-1)/factorial(i-1);            
    end
end
 
%Al
% b = zeros(7,1);
% b(2) = 1;

cl = Al\b;

cr = -cl;
cr = cr(end:-1:1);

dUdt = zeros(m,n);
for j = 1:n
    if(j > 4 && j < n-3)
    dUdt(:,j)= U(:,j-4:j+4)*cc/k;
    elseif(j<=4)
        [j n];
        dUdt(:,j) = U(:,j:j+8)*cl/k;
    elseif(j>=(n-3))
       
        dUdt(:,j) = U(:,j-8:j)*cr/k;
    end
    
end
% U
% dUdt
% t = 0:0.1:2;
% plot(t,dUdt,t,-sin(t),t,cos(t))




end

