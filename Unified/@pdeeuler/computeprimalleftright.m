function computeprimalleftright(obj  )
%COMPUTEPRIMALLEFTRIGHT Summary of this function goes here
%   Detailed explanation goes here
Z = obj.convVreconp;
N = obj.nCells;
h = obj.cellWidths;


if(obj.NLfluxtype==4)
   
U = obj.convSoln;

U1r = zeros(N+2,1);
U1l = zeros(N+2,1);
U2r = zeros(N+2,1);
U2l = zeros(N+2,1);
U3r = zeros(N+2,1);
U3l = zeros(N+2,1);

p = obj.pOrder;
obj.pOrder = obj.qOrder;
obj.computeprimalpseudo();

order  = obj.pOrder;
 Z = obj.unstructuredrecon(U,order,'solution');
for i = 2:N+1
    for k = 1:order
    U1r(i) = U1r(i)+ Z(k,i)*(h(i)/2)^(k-1);
    U1l(i) = U1l(i)+ Z(k,i)*(-h(i)/2)^(k-1);
    U2r(i)   = U2r(i)+ Z(k+order,i)*(h(i)/2)^(k-1);
    U2l(i)   = U2l(i)+ Z(k+order,i)*(-h(i)/2)^(k-1);
    U3r(i)   = U3r(i)+ Z(k+2*order,i)*(h(i)/2)^(k-1);
    U3l(i)   = U3l(i)+ Z(k+2*order,i)*(-h(i)/2)^(k-1);
    end
end
obj.convUleft = [U1l U2l U3l];
obj.convUright = [U1r U2r U3r];
 
    

obj.pOrder = p;
obj.computeprimalpseudo();

   return; 
end


% error('1')
p = obj.pOrder;
obj.pOrder = obj.qOrder;
obj.computeprimalpseudo();
%

primrhor = zeros(N+2,1);
primrhol = zeros(N+2,1);
primur = zeros(N+2,1);
primul = zeros(N+2,1);
primPr = zeros(N+2,1);
primPl = zeros(N+2,1);
order  = obj.pOrder;
Z = obj.unstructuredrecon(obj.convSolutionV,order,'solution');
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

%
obj.pOrder = p;
obj.computeprimalpseudo();



end