clearvars -except errerr10 errerr20 errerr40 errerr10b errerr20b errerr40b berrerr2 berrerr4 berrerr8 berrerr berrerr4 berrerr8 berrerrnon2 berrerrnon4 berrerrlin2 berrerrlin4 err2b cverr2b err4b x exacterr ee
close all
A = [2 2 4;2 2 6;2 4 4;2 4 6; 2 6 4; 2 6 6; 4 2 6; 4 4 6; 4 6 6];
B = ['2 2 4';'2 2 6';'2 4 4';'2 4 6'; '2 6 4'; '2 6 6'; '4 2 6'; '4 4 6'; '4 6 6'];

A = [2 2 5; 2 3 4; 2 3 5; 3 3 5; 3 4 5; 3 5 5 ; 4 2 5; 4 4 5; 4 5 6];
A = [3 2 4; 3 2 5; 3 2 6; 3 3 4; 4 3 5; 4 3 6; 4 2 5; 4 2 6; 5 2 6];
A = [2 2 3; 2 2 4; 2 2 5; 2 2 6; 2 3 3; 2 3 4; 2 3 5; 2 3 6; 2 4 3; 2 4 4; 2 4 5; 2 4 6; 2 5 3; 2 5 4; 2 5 5; 2 5 6; 2 6 3; 2 6 4; 2 6 5; 2 6 6; 3 2 4; 3 2 5; 3 2 6; 3 3 4; 3 3 5; 3 3 6; 3 4 4; 3 4 5; 3 4 6; 3 5 4; 3 5 5; 3 5 6; 3 6 4; 3 6 5; 3 6 6; 4 2 5; 4 2 6; 4 3 5; 4 3 6; 4 4 5; 4 4 6; 4 5 5 ; 4 5 6; 4 6 5; 4 6 6; 5 2 6; 5 3 6 ;5 4 6 ; 5 5 6; 5 6 6  ]
B = ['(2,2,3)'; '(2,2,4)'; '(2,2,5)'; '(2,2,6)'; '(2,3,3)'; '(2,3,4)'; '(2,3,5)'; '(2,3,6)'; '(2,4,3)'; '(2,4,4)'; '(2,4,5)'; '(2,4,6)'; '(2,5,3)'; '(2,5,4)'; '(2,5,5)'; '(2,5,6)'; '(2,6,3)'; '(2,6,4)'; '(2,6,5)'; '(2,6,6)';];% '3 2 4'; '3 2 5'; '3 2 6'; '3 3 4'; '3 3 5'; '3 3 6'; '3 4 4'; '3 4 5'; '3 4 6'; '3 5 4'; '3 5 5'; '3 5 6'; '3 6 4'; '3 6 5'; '3 6 6'; '4 2 5'; '4 2 6'; '4 3 5'; '4 3 6'; '4 4 5'; '4 4 6'; '4 5 5' ; '4 5 6'; '4 6 5'; '4 6 6'; '5 2 6'; '5 3 6' ;'5 4 6' ; '5 5 6'; '5 6 6'  ]

[m,n] = size(A);
N = 40
% figure
% hold on
% for j = 1:m
%    [err4bb(j),x,cverr4b(j),exacterr(:,j),ee(:,j)] = errordriver(N,A(j,1),A(j,2),A(j,3),1/3,100,.3,7,'Advection');
%  
%   %[errerr1b(j),x,cverr,exacterr(:,j),ee(:,j)] = errordriver(N,A(j,1),A(j,2),A(j,3),1,100,10,7,'Poisson');
%  %%% [errerr4(j),x,cverr,exacterr(:,j),ee(:,j)] = errordriver(N,A(j,1),A(j,2),A(j,3),0,100,.1,7,'Poisson');
%  
% %  [berrerrlin8(j),x,cverr,exacterr(:,j),ee(:,j)] = errordriver(N,A(j,1),A(j,2),A(j,3),0,100,.1,7,'Burgers');
% close all
% 
% end

load('errorarrayplot.mat');

%    [err4bb(20),x,cverr4b(20),exacterr(:,20),ee(:,20)] = errordriver(N,2,6,6,1/3,100,.3,7,'Advection');

h=figure
set(gca,'FontSize',14)
for k = 1:20
subplot(5,4,k)
plot(x,exacterr(:,k),'o',x,ee(:,k),'*');
ylim([-1e-2 1e-2])
xlabel(B(k,:),'FontSize',14)
% ylabel('solution error')
   set(gca,'FontSize',14)
if (k==1)
   legend('exact error','computed error'); 

   
   
end

hold on
% subplot(2,1,k)
% plot(x,exacterr(:,k),'o',x,ee(:,2),'*');


   if(k==9)
     ylabel('$$\mathbf{\mathcal{I}^h\epsilon},\bf{\epsilon}_h$$','Interpreter','Latex','FontSize',24) 
   end


end
% annotation(h,'textbox',...
%     [0.0504791666666667 0.532837474674436 0.0322916666666667 0.0351288056206089],...
%     'Interpreter','latex',...
%     'String',{'$$\mathbf{\mathcal{I}^h\epsilon},\bf{\epsilon}_h$$'},...
%     'FontWeight','bold',...
%     'FontSize',24,...
%     'FontName','AlArabiya',...
%     'EdgeColor',[1 1 1]);

annotation(h,'textbox',...
    [0.506208333333333 0.0348869634744324 0.0296875 0.0304449648711944],...
    'Interpreter','latex',...
    'String',{'$$x$$'},...
    'FontWeight','bold',...
    'FontSize',28,...
    'FontName','AlArabiya',...
    'EdgeColor',[1 1 1]);






type = {'xb', '^k','vr','sm','dg'}
tab = [ 3 2 2 2 3 2 2 2 3 4 2 2 3 2 5 2 3 2 2 6 ]';
g=figure
kk=1;
for k = 1:20
subplot(5,4,k)
v = plot(x,exacterr(:,k)-ee(:,k),type{tab(k)-1},'LineWidth',2);



 if(~((mod(k,4)==1)|| (k>=10 && mod(k,5)==0)))
   ylim([-1e-3 1e-3]) 
   set(gca,'YTick',[-1e-3 1e-3])
 else
     
     yt = get(gca,'ytick');
     set(gca,'YTick',[yt(1) yt(3)])
 end
 set(gca,'XTick',[0 1])

if(k==2 || k== 5 || k==10 || k==15 || k==20)
lin(kk) = v;
kk = kk+1;

end
   

    set(gca,'FontSize',24)
 
% % % xlabel(B(k,:),'FontSize',14)


% if (k==1)
%    legend('exact error','computed error'); 
% end
hold on
% subplot(2,1,k)
% plot(x,exacterr(:,k),'o',x,ee(:,2),'*');
   if(k==9)
     ylabel('$$\mathbf{\mathcal{I}^h\epsilon-{\epsilon}_h}$$','Interpreter','Latex','FontSize',35) 
   end

end

subplot(5,4,20)
v= plot(x,exacterr(:,20)-ee(:,20),'Color',[153/255, 76/255, 0],'LineStyle', 'd','LineWidth',2);
lin(kk-1) = v;
% annotation(g,'textbox',...
%     [0.0504791666666667 0.532837474674436 0.0322916666666667 0.0351288056206089],...
%     'Interpreter','latex',...
%     'String',{'$$\mathbf{\mathcal{I}^h\epsilon-{\epsilon}_h}$$'},...
%     'FontWeight','bold',...
%     'FontSize',24,...
%     'FontName','AlArabiya',...
%     'EdgeColor',[1 1 1]);

annotation(g,'textbox',...
    [0.506208333333333 0.0348869634744324 0.0296875 0.0304449648711944],...
    'Interpreter','latex',...
    'String',{'$$x$$'},...
    'FontWeight','bold',...
    'FontSize',35,...
    'FontName','AlArabiya',...
    'EdgeColor',[1 1 1]);

hL = legend([lin],{'$$\mathcal{O}(h^\mathbf{2})$$','$$\mathcal{O}(h^\mathbf{3})$$','$$\mathcal{O}(h^\mathbf{4})$$','$$\mathcal{O}(h^\mathbf{5})$$','$$\mathcal{O}(h^\mathbf{6})$$'},'Interpreter','Latex','FontSize',30,'Orientation','Horizontal');



%zoom in 
fig = figure

F = 21;
g = 3;

       f(1)= subplot(2,2,1);
linn(1)=plot(x,exacterr(:,5)-ee(:,5),'^','Color','Black','LineWidth',g);
         xlabel('(2,3,3)','FontSize',F+2);
                set(f(1), 'fontsize', F);
    set(findobj(f(1),'Type','text'),'FontSize',  F);
                 
       f(2)= subplot(2,2,2);
linn(2)=plot(x,exacterr(:,10)-ee(:,10),'v','Color','Red','LineWidth',g);
         xlabel('(2,4,4)','FontSize',F+2);
                set(f(2), 'fontsize', F);
    set(findobj(f(2),'Type','text'),'FontSize',  F);
                 
           f(3)= subplot(2,2,3);
           
           
           
linn(3)=plot(x,exacterr(:,15)-ee(:,15),'s','Color','Magenta','LineWidth',g);
         xlabel('(2,5,5)','FontSize',F+2);
                set(f(3), 'fontsize', F);
    set(findobj(f(3),'Type','text'),'FontSize',  F);
%              ylabel('\hspace*{2.5cm}$$\mathbf{\epsilon_p-{\epsilon}_{pq}}$$','Interpreter','Latex','FontSize',30)       
           f(4)= subplot(2,2,4);
linn(4)=plot(x,exacterr(:,20)-ee(:,20),'d','Color',[153/255, 76/255, 0],'LineWidth',g);
         xlabel('(2,6,6)','FontSize',F+2);
                set(f(4), 'fontsize', F);
    set(findobj(f(4),'Type','text'),'FontSize',  F);
                 
gL = legend([linn],{'$$\mathcal{O}(h^\mathbf{3})$$','$$\mathcal{O}(h^\mathbf{4})$$','$$\mathcal{O}(h^\mathbf{5})$$','$$\mathcal{O}(h^\mathbf{6})$$'},'Interpreter','Latex','FontSize',22,'Orientation','Horizontal');

annotation(fig,'textbox',...
    [0.506208333333333 0.0348869634744324 0.0296875 0.0304449648711944],...
    'Interpreter','latex',...
    'String',{'$$x$$'},...
    'FontWeight','bold',...
    'FontSize',25,...
    'FontName','AlArabiya',...
    'EdgeColor',[1 1 1]);




% for k = 21:35
% subplot(5,3,k-20)
% plot(x,exacterr(:,k),'o',x,ee(:,k),'*');
% 
% xlabel(B(k,:))
% if (k==1)
%    legend('exact error','computed error'); 
% end
% hold on
% % subplot(2,1,k)
% % plot(x,exacterr(:,k),'o',x,ee(:,2),'*');
% end
% 
% figure
% for k = 21:35
% subplot(5,3,k-20)
% plot(x,exacterr(:,k)-ee(:,k),'x');
% 
% xlabel(B(k,:))
% % if (k==1)
% %    legend('exact error','computed error'); 
% % end
% hold on
% % subplot(2,1,k)
% % plot(x,exacterr(:,k),'o',x,ee(:,2),'*');
% end
% 
