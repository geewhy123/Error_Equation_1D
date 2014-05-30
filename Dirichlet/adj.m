clear all
close all
load('adj.mat')
load('dual.mat')


N = 10;
h = (1/N)*ones(N+2,1);
h(1) = NaN;
h(N+2) = NaN;
p=2;
global AD
AD = computepseudo(N,x,h,p);   
[Z] = unstructuredrecon(u,x,h,N,0,0,p)
 [er]=reconplot(x,h,N,p,Z);