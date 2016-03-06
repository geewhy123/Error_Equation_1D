function [ R ] = computeres(obj,u,time,r)%u,x,h,N,f,r,phys,time,gsp )
%COMPUTERES Summary of this function goes here
%   Detailed explanation goes here
global TEND
if(abs(time-round(time))<1e-10)
    time = round(time);
elseif(abs(time-TEND) < 1e-10)
    time = TEND;
end

N = obj.nCells;
phys = obj.physics;

eqn = 'residual';
if(r==obj.pOrder)
    eqn = 'solution';
end
Z = obj.unstructuredrecon(u,r,eqn);

R= obj.computefluxintegral(Z,'residual');
% Z
% R
% error('1')
for i = 2:N+1       
    global dUdt
    global KK
    ut = dUdt(i,round(time/KK)+1);
    
    Ut(i) = ut;
    
    
    % % % % R(i)= -ut+(right(i)-left(i))/h(i)-f(i);
    R(i) = R(i) - ut;
    
    
    if(strcmp(phys,'Burgers')==1)
        error('1')
        R(i) = -ut-(right(i)^2-left(i)^2)/(2*h(i))-f(i);
    end
    
    
end
