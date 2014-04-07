function [ upr,upl,phi]  = timestep( equation, Z,f,k,h,N,p,phys,time,Rsp,Zu,val)
%   Detailed explanation goes here

switch equation
    case 'solution'
        [upr,upl,phi]=reconfluxsoln( Z,f,h,N,p,phys,time);
    case 'error'
        [upr,upl,phi]=reconfluxerror( Z,f,k,h,N,p,phys,time,Rsp,Zu,val);
    otherwise
        error('2')
    

end

