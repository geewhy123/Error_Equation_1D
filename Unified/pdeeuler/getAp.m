function [ Ap ] = getAp( xx )
%GETAP Summary of this function goes here
%   Detailed explanation goes here

Ae = 0.4;
At = 0.2;
Ap =  0;%(50/9)*(Ae-At)*(xx-2/5); 
end

