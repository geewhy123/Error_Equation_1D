close all
clear all
% % load('tauN.mat')
% % 
% % tau3264 = tau32;
% % N = 32;
% % for j = 2:N+1
% %    h12 = h64(2*(j-1))
% %    h23 = h64(2*(j-1)+1);
% %    tau3264(j) = (h12*tau64(2*(j-1))+h23*tau64(2*(j-1)+1))/(h12+h23); 
% % end
% % plot(x32,tau32,'-+',x64,tau64,'-*',x32,tau3264,'o')
% % set(gca,'FontSize',15)
% % legend('N=32','N=64','N=64 interpolated back to N=32 by weighted average')
% % xlabel('x')
% % ylabel('\tau')
% % 
% % save('tauN.mat','tau3264','-append')

% MG

load('te.mat')


u2040b = u20b;

for j = 2:length(u2040b)-1
    k = 2*j-2;
   
   u2040b(j) = (h40b(k)*u40b(k)+h40b(k+1)*u40b(k+1))/(h40b(k)+h40b(k+1));
end


subplot(211)
plot(x20b,u20b,'o',x20b,u2040b,'*')
legend('actual coarse soln','coarse soln based on MG fine solution')
subplot(212)
plot(x20b,tau20b,'o',x20b,tau2040b,'*')
legend('actual coarse t.e.','t.e. based on MG fine solution')
% plot(x20,u2040,'*',x20,u20,'o')
save('te.mat','u2040b','-append')
% save('te.mat','u2040','-append')