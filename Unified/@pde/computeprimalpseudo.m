function  computeprimalpseudo(obj)
%COMPUTEPSEUDO Summary of this function goes here
%   Detailed explanation goes here
% 
    p = obj.pOrder;
    switch p
    case 6
        fprintf('6th order need more testing')
%               assert(0)
        [AD,AA]=obj.computepseudo(p);
    otherwise
        [AD,AA]=obj.computepseudo(p);
    end

    obj.primalPI = AD;
    obj.primalRM = AA;
end








