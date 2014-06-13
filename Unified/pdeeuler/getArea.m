function [ A] = getArea( xx )
%GETAREA Summary of this function goes here
%   Detailed explanation goes here

Ae = 0.4;
At = 0.2;
A =  1;%(25/9)*(Ae-At)*(xx-2/5).^2+At; 

end

