% clear all
close all
p = [2 ; 3; 4; 5];
q = [3 ; 4; 5; 6];
r = q;
t = [2 ;3 ;4 ;7];
m = 4;

N = [20 40 80];
figure
hold on


for i = 1:m
    for j = 1:m
        for k = 1:m
            for l = 1:length(N)
                pp = p(i);
                qq = q(j);
                rr = r(j);
                tt = t(k);
                ee(i,j,k,l) = errordriver(N(l),pp,qq,rr,1/3,'P',0,'P',0,1,tt,'Advection','TimeAccurate');
            end
        end
    end
end

save('unsteadytestfull.mat','ee')
% load('unsteadytest.mat')
str = ['-*' ;'-o' ;'-^'];
for i = 1:5
subplot(2,3,i)
for j = 1:3
loglog(N.^-1,ee(j,:,i),str(j,:))
log(ee(j,end,i)/ee(j,end-1,i))/log(2)
hold on
end
grid on
end