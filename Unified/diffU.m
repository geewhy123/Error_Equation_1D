    function [ dUdt] = diffU( U ,k,klast)
%DIFFU Summary of this function goes here
%   Detailed explanation goes here

[m,n] = size(U);

% k
% size(U)
% U(:,end-2:end)
% error('1')

% 6th order
%A = [1 1 1 1 1 1 1 1 1 ; -4 -3 -2 -1 0 1 2 3 4; 
%b = [0 1 0]';

%centered
nt = 8;

for j = 1:nt+1
    for i = 1:nt+1
     A(i,j) = (j-(nt/2+1))^(i-1)/factorial(i-1);            
    end
end
                
    b = zeros(nt+1,1);
    b(2) = 1;
cc = A\b;

for j = 1:nt+1
    for i = 1:nt+1
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
    if(j > nt/2 && j < n-(nt/2-1))
    dUdt(:,j)= U(:,j-nt/2:j+nt/2)*cc/k;
    elseif(j<=nt/2)
        [j n];
        dUdt(:,j) = U(:,j:j+nt)*cl/k;
    elseif(j>=(n-(nt/2-1)))
       
        dUdt(:,j) = U(:,j-nt:j)*cr/k;
    end
    
end
% U
% dUdt
% t = 0:0.1:2;
% plot(t,dUdt,t,-sin(t),t,cos(t))




end

