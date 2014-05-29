function  computerespseudo(obj)
%COMPUTEPSEUDO Summary of this function goes here
%   Detailed explanation goes here
% 
r = obj.rOrder;
switch r
           case 6
               fprintf('6th order doesnt work yet')
%                assert(0)
        [AD,AA]=obj.computepseudo(r);
    otherwise
        [AD, AA]=obj.computepseudo(r);
        
        
end



obj.resPI = AD;
obj.resRM = AA;

end









