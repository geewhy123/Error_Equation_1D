function [ R ,uxx,Z] = computeres(u,x,h,N,f )
%COMPUTERES Summary of this function goes here
%   Detailed explanation goes here

[error,Z]=unstructuredrecon4(u,x,h,N,1,exp(-1));
R = zeros(N+2,1);
uxx = zeros(N+2,1);
Z
for i = 2:N+1
    hi = h(i);
    xi = x(i);
   x2bi = (1/hi)*( ((xi+hi/2)-xi)^3/3-((xi-hi/2)-xi)^3/3);
   
    uxx(i) = 2*Z(3,i)+12*Z(5,i)*x2bi;
   
  %uxx(i) = 2*Z(3,i);
%   uxx(i) = 2*Z(3,i);%+6*Z(4,i)*(x-
   R(i) = uxx(i) - f(i);
end
uxx(N/2)
end

