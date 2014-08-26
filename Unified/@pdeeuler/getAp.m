function [ Ap ] = getAp( obj,xx )
%GETAP Summary of this function goes here
%   Detailed explanation goes here

Ae = 0.4;
At = 0.2;

if(obj.areatype == 0)
Ap =  0;
elseif(obj.areatype == 1)
Ap = (50/9)*(Ae-At)*(xx-2/5); 
elseif(obj.areatype ==2)
    a1 = 0.1;
    Ap = -2*pi*a1*sin(2*pi*xx);
    end
end

