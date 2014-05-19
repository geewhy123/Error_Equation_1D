function [error, Z] = unstructuredrecon1( u,x,h,N,BCLeft,u0,BCRight,u1 )
%UNSTRUCTUREDRECON3 Summary of this function goes here
%   Detailed explanation goes here
Z = zeros(2,N+2);
error = 0;

%4th order recon

if(BCLeft == 'P' && BCRight == 'P')
for i = 2:N+1


switch i
    case 2
 Y= recon1(x(i),h(i),u(i),x(N+1),h(N+1),u(N+1),x(i+1),h(i+1),u(i+1),x(N),h(N),u(N),x(i+2),h(i+2),u(i+2),i);
    case 3
  Y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(N+1),h(N+1),u(N+1),x(i+2),h(i+2),u(i+2),i);
    case N
 Y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(2),h(2),u(2),i);
    case N+1
   Y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(2),h(2),u(2),x(i-2),h(i-2),u(i-2),x(3),h(3),u(3),i);
    otherwise
  Y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),i);
 
end


y=Y;
 Z(:,i) = y;



end

elseif(BCLeft=='D' && BCRight == 'D')
    for i = 2:N+1

switch i
    case 2
         Y= reconboundary1(x(i),h(i),u(i),x(i+1),h(i+1),u(i+1),x(i+2),h(i+2),u(i+2),x(i+3),h(i+3),u(i+3),x(i+4),h(i+4),u(i+4),u0,'left',i,NaN);
%  Y= recon1(x(i),h(i),u(i),x(N+1),h(N+1),u(N+1),x(i+1),h(i+1),u(i+1),x(N),h(N),u(N),x(i+2),h(i+2),u(i+2),i);
    case 3
         Y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i+2),h(i+2),u(i+2),x(i+3),h(i+3),u(i+3),i);
%   Y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(N+1),h(N+1),u(N+1),x(i+2),h(i+2),u(i+2),i);
    case N
        Y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i-3),h(i-3),u(i-3),i);
%  Y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(2),h(2),u(2),i);
    case N+1
        Y= reconboundary1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i-2),h(i-2),u(i-2),x(i-3),h(i-3),u(i-3),x(i-4),h(i-4),u(i-4),u1,'right',i,NaN);
%    Y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(2),h(2),u(2),x(i-2),h(i-2),u(i-2),x(3),h(3),u(3),i);
    otherwise
      Y= recon1(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),i);
 

  
end
y=Y;
 Z(:,i) = y;



end

    
end

