function  computeerrorpseudo(obj)
%COMPUTEPSEUDO Summary of this function goes here
%   Detailed explanation goes here
% 
    p = obj.qOrder;
    switch p
    case 6
        fprintf('6th order doesnt work yet')
%                assert(0)
        [AD,AA]=obj.computepseudo(p);
    otherwise
        [AD,AA]=obj.computepseudo(p);
    end

    obj.errorPI = AD;
    obj.errorRM = AA;
end







