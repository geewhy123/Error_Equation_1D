function [ res ] = getRes( time ,k,i,sp)
%GETRES Summary of this function goes here
%   Detailed explanation goes here


global R
global M
nSteps = M;

%nSteps = 50;
T=(0:1:nSteps)*(k);


%sp = spapi(6,T,R(i,:));




res=fnval(sp,time);

end

