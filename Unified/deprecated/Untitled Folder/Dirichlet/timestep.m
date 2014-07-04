function [ phi]  = timestep( equation, Z,f,k,h,N,p,phys,time,Rsp,Zu,val)
%   Detailed explanation goes here

switch equation
    case 'solution'
        phi=reconfluxsoln( Z,f,h,N,p,phys,time);
    case 'error'
        phi=reconfluxerror( Z,f,k,h,N,p,phys,time,Rsp,Zu,val);
    otherwise
        error('2')
    

end

