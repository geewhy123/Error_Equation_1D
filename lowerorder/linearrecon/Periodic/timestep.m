function [ upr,upl,phi]  = timestep( equation, u,Z,f,k,h,N,p,phys,uder,j,time,gsp,Rsp,Zu)
%   Detailed explanation goes here

switch equation
    case 'solution'
        [upr,upl,phi]=reconfluxsoln( u,Z,f,k,h,N,p,phys,uder,j,time,gsp,Rsp,Zu);
    case 'error'
        [upr,upl,phi]=reconfluxerror( u,Z,f,k,h,N,p,phys,uder,j,time,gsp,Rsp,Zu);
    otherwise
        error('2')
    

end

