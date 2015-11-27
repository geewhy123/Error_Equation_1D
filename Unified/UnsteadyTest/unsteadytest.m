% clear all

%burgerstest

% error('1')
close all
p = [2 ; 3; 4; 5];
q = [3 ; 4; 5; 6];
r = q;
t = [2 ;3 ;4 ;7];
t = ['irk1';'irk2';'irk4'];
m = 4;

N = [20 40 80];
figure
hold on

% load('unsteadytestfullprofilesdiffusionimplicit.mat')
% for i = 1:m
%     for j = 1:m
%         for k = 1:3
%             for l = 1:length(N)
%                 pp = p(i);
%                 qq = q(j);
%                 rr = r(j);
%                 tt = t(k,:);
% %                 if(pp ~= 5 && qq ~= 5)
% %                    continue 
% %                 end
%                 [pp qq rr ]
%                 tt
% %                 ee(i,j,k,l) = errordriver(N(l),pp,qq,rr,1/3,'P',0,'P',0,1,tt,'Advection','TimeAccurate');
% if(l == 1)
% [errerr(i,j,k,l),x,cverr(i,j,k,l),exacterr(i,j,k,:),ee(i,j,k,:)] = errordriver(N(l),pp,qq,rr,1/3,'D',0,'D',0,0.1,tt,'Poisson','TimeAccurate');
% % [errerr(i,j),x,cverr(i,j),exacterr(i,j,:),ee(i,j,:)] = errordriver(N(l),pp,qq,rr,1/3,'P',0,'P',0,0.1,tt,'Poisson','TimeAccurate');
% else
%     errerr(i,j,k,l) = errordriver(N(l),pp,qq,rr,1/3,'D',0,'D',0,0.1,tt,'Poisson','TimeAccurate');
% end
%             end
%         end
%     end
% end

% save('unsteadytestfull.mat','ee')
% save('unsteadytestfullprofilesdiffusionimplicitdirichlet.mat','errerr','x','cverr','exacterr','ee')
% error('1')
% load('unsteadytest.mat')
% load('unsteadytestfull.mat')
% load('unsteadytestfullprofiles.mat')
% load('unsteadytestfullprofilesdiffusion3.mat')
% load('unsteadytestfullprofilesdiffusionimplicit.mat')
load('unsteadytestfullprofilesdiffusionimplicitdirichlet.mat')
str = ['-*' ;'-o' ;'-^';'-v'];
% for i = 1:5
% subplot(2,3,i)
% for j = 1:3
% loglog(N.^-1,ee(j,:,i),str(j,:))
% log(ee(j,end,i)/ee(j,end-1,i))/log(2)
% hold on
% end
% grid on
% end

T = 3;

figure
for i = 1:length(p)
    for j = 1:length(q)
       subplot(length(p),length(q),4*(i-1)+j)
       for k = 1:length(N)
          data(k) = errerr(i,j,T,k); 
       end
       loglog(N.^-1,data,str(j,:))
       order(i,j) = log(data(2)/data(3))/log(2);
    end 
end
order


figure

for i = 1:length(p)
    for j = 1:length(q)
        for m = 1:length(x);
           err(m) = exacterr(i,j,T,m)-ee(i,j,T,m); 
        end
       subplot(length(p),length(q),4*(i-1)+j)
       plot(x,err,'-*')
    end 
end

