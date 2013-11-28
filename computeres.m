function [ R ,uxx,Z] = computeres(u,x,h,N,f )
%COMPUTERES Summary of this function goes here
%   Detailed explanation goes here

[error,Z]=unstructuredrecon3(u,x,h,N,1,exp(-1));
R = zeros(N+2,1);
uxx = zeros(N+2,1);

for i = 2:N+1
   uxx(i) = 2*Z(3,i);%+6*Z(4,i)*(x-
   R(i) = uxx(i) - f(i);
end

end

