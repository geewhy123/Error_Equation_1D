function [ AD ] = computepseudo(N,x,h,p)
%COMPUTEPSEUDO Summary of this function goes here
%   Detailed explanation goes here
% 
switch p
    case 2
         [AD]=computepseudo1(N,x,h);
    case 3
         [AD]=computepseudo2(N,x,h);
    case 4
        [AD]=computepseudo3(N,x,h);
    case 5
        [AD]=computepseudo4(N,x,h);
           case 6
        [AD]=computepseudo5(N,x,h);
    otherwise
        error('1')
end





end

