function [ res ] = getRes( time ,k,i,sp)
%GETRES Summary of this function goes here
%   Detailed explanation goes here


global R
global M
nSteps = M;

%nSteps = 50;
T=(0:1:nSteps)*(k);


%sp = spapi(6,T,R(i,:));



% tic
res=fnval(sp,time);
% toc
% error('1')
%res
%error('1')
end

