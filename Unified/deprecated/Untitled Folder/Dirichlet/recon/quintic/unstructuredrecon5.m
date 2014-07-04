function [error, Z] = unstructuredrecon5( u,x,h,N,u0,u1 )
%UNSTRUCTUREDRECON3 Summary of this function goes here
%   Detailed explanation goes here
Z = zeros(6,N+2);
error = 0;


for i = 2:N+1
switch i
    case 2
         Y= reconboundary5(x(i),h(i),u(i),x(i+1),h(i+1),u(i+1),x(i+2),h(i+2),u(i+2),x(i+3),h(i+3),u(i+3),x(i+4),h(i+4),u(i+4),x(i+5),h(i+5),u(i+5),x(i+6),h(i+6),u(i+6),u0,'left',i);
%  Y= recon5(x(i),h(i),u(i),x(N+1),h(N+1),u(N+1),x(i+1),h(i+1),u(i+1),x(N),h(N),u(N),x(i+2),h(i+2),u(i+2),x(N-1),h(N-1),u(N-1),x(i+3),h(i+3),u(i+3),i);
    case 3
          Y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i+2),h(i+2),u(i+2),x(i+3),h(i+3),u(i+3),x(i+4),h(i+4),u(i+4),x(i+5),h(i+5),u(i+5),i);
%   Y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(N+1),h(N+1),u(N+1),x(i+2),h(i+2),u(i+2),x(N),h(N),u(N),x(i+3),h(i+3),u(i+3),i);
    case 4
          Y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),x(i+3),h(i+3),u(i+3),x(i+4),h(i+4),u(i+4),i);
%          Y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),x(N+1),h(N+1),u(N+1),x(i+3),h(i+3),u(i+3),i);
    case N-1
         Y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),x(i-3),h(i-3),u(i-3),x(i-4),h(i-4),u(i-4),i);
%          Y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),x(i-3),h(i-3),u(i-3),x(2),h(2),u(2),i);
    case N
           Y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i-3),h(i-3),u(i-3),x(i-4),h(i-4),u(i-4),x(i-5),h(i-5),u(i-5),i);
%  Y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(2),h(2),u(2),x(i-3),h(i-3),u(i-3),x(3),h(3),u(3),i);
    case N+1
          Y= reconboundary5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i-2),h(i-2),u(i-2),x(i-3),h(i-3),u(i-3),x(i-4),h(i-4),u(i-4),x(i-5),h(i-5),u(i-5),x(i-6),h(i-6),u(i-6),u1,'right',i);
%    Y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(2),h(2),u(2),x(i-2),h(i-2),u(i-2),x(3),h(3),u(3),x(i-3),h(i-3),u(i-3),x(4),h(4),u(4),i);
    otherwise
  Y= recon5(x(i),h(i),u(i),x(i-1),h(i-1),u(i-1),x(i+1),h(i+1),u(i+1),x(i-2),h(i-2),u(i-2),x(i+2),h(i+2),u(i+2),x(i-3),h(i-3),u(i-3),x(i+3),h(i+3),u(i+3),i);
 
   

  
end

 y=Y;
   Z(:,i) = y;


end

