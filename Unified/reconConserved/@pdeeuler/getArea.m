function [ A] = getArea(obj, xx )
%GETAREA Summary of this function goes here
%   Detailed explanation goes here

Ae = 0.4;
At = 0.2;

if(obj.areatype == 0)
A =  1;
elseif(obj.areatype == 1)
A = (25/9)*(Ae-At)*(xx-2/5).^2+At; 
elseif(obj.areatype == 2)
    a0 = 0.3;
    a1 = 0.1;
    A = a0+a1*cos(2*pi*xx);
        
end

end

