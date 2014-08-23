function computeprimalleftright(obj  )
%COMPUTEPRIMALLEFTRIGHT Summary of this function goes here
%   Detailed explanation goes here

Z = obj.convVreconp;
N = obj.nCells;
h = obj.cellWidths;

primrhor = zeros(N+2,1);
primrhol = zeros(N+2,1);
primur = zeros(N+2,1);
primul = zeros(N+2,1);
primPr = zeros(N+2,1);
primPl = zeros(N+2,1);
order  = obj.pOrder;
for i = 2:N+1
    for k = 1:order
    primrhor(i) = primrhor(i)+ Z(k,i)*(h(i)/2)^(k-1);
    primrhol(i) = primrhol(i)+ Z(k,i)*(-h(i)/2)^(k-1);
    primur(i)   = primur(i)+ Z(k+order,i)*(h(i)/2)^(k-1);
    primul(i)   = primul(i)+ Z(k+order,i)*(-h(i)/2)^(k-1);
    primPr(i)   = primPr(i)+ Z(k+2*order,i)*(h(i)/2)^(k-1);
    primPl(i)   = primPl(i)+ Z(k+2*order,i)*(-h(i)/2)^(k-1);
    end
end



obj.convVleft = [primrhol primul primPl];
obj.convVright = [primrhor primur primPr];
obj.convUleft = zeros(N+2,3);
for i = 2:N+1
[obj.convUleft(i,1),obj.convUleft(i,2),obj.convUleft(i,3)] = toconservedvars(primrhol(i),primul(i),primPl(i));
end
obj.convUright = zeros(N+2,3);
for i = 2:N+1
[obj.convUright(i,1),obj.convUright(i,2),obj.convUright(i,3)] = toconservedvars(primrhor(i),primur(i),primPr(i));
end


end