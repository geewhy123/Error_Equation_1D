close all
clear all

% % % single method comparison
% N = 20;
% figure
% 
% [errerr2,x1,cverr2,exacterr1,ee1,te] = errordriver(N,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',1)
% [errerr2,x2,cverr2,exacterr2,ee2,te] = errordriver(N,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',2)
% [errerr2,x3,cverr2,exacterr3,ee3,te] = errordriver(N,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',3)
% [errerr2,x4,cverr2,exacterr4,ee4,te] = errordriver(N,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',4)
% 
% close all
% figure 
% set(gca,'FontSize',15)
% subplot(3,4,1)
% plot(x1,ee1(:,1),'*',x1,exacterr1(:,1),'o','LineWidth',2);
% title('Method 1','FontSize',20)
% subplot(3,4,5)
% plot(x1,ee1(:,2),'*',x1,exacterr1(:,2),'o','LineWidth',2);
% ylabel('Error in Primitive','FontSize',18)
% subplot(3,4,9)
% plot(x1,ee1(:,3),'*',x1,exacterr1(:,3),'o','LineWidth',2);
% xlabel('x','FontSize',18)
% 
% subplot(3,4,2)
% plot(x2,ee2(:,1),'*',x2,exacterr2(:,1),'o','LineWidth',2);
% title('Method 2','FontSize',20)
% subplot(3,4,6)
% plot(x2,ee2(:,2),'*',x2,exacterr2(:,2),'o','LineWidth',2);
% ylabel('Error in Conserved','FontSize',18)
% subplot(3,4,10)
% plot(x2,ee2(:,3),'*',x2,exacterr2(:,3),'o','LineWidth',2);
% xlabel('x','FontSize',18)
% 
% subplot(3,4,3)
% plot(x3,ee3(:,1),'*',x3,exacterr3(:,1),'o','LineWidth',2);
% title('Method 3','FontSize',20)
% subplot(3,4,7)
% plot(x3,ee3(:,2),'*',x3,exacterr3(:,2),'o','LineWidth',2);
% ylabel('Error in Conserved','FontSize',18)
% subplot(3,4,11)
% plot(x3,ee3(:,3),'*',x3,exacterr3(:,3),'o','LineWidth',2);
% xlabel('x','FontSize',18)
% 
% subplot(3,4,4)
% plot(x4,ee4(:,1),'*',x4,exacterr4(:,1),'o','LineWidth',2);
% title('Method 4','FontSize',20)
% subplot(3,4,8)
% plot(x4,ee4(:,2),'*',x4,exacterr4(:,2),'o','LineWidth',2);
% ylabel('Error in Conserved','FontSize',18)
% subplot(3,4,12)
% plot(x4,ee4(:,3),'*',x4,exacterr4(:,3),'o','LineWidth',2);
% xlabel('x','FontSize',18)
% h = legend('Error Estimate','Exact Error');
% set(h,'FontSize',13)
% error('1')

% % % BW plot


N = 20;
%   Q = [2 2 3; 2 2 4; 2 2 5; 2 2 6; 2 3 3; 2 3 4; 2 3 5; 2 3 6; 2 4 3; 2 4 4; 2 4 5; 2 4 6; 2 5 3; 2 5 4; 2 5 5; 2 5 6; 2 6 3; 2 6 4; 2 6 5; 2 6 6; 3 2 4; 3 2 5; 3 2 6; 3 3 4; 3 3 5; 3 3 6; 3 4 4; 3 4 5; 3 4 6; 3 5 4; 3 5 5; 3 5 6; 3 6 4; 3 6 5; 3 6 6; 4 2 5; 4 2 6; 4 3 5; 4 3 6; 4 4 5; 4 4 6; 4 5 5 ; 4 5 6; 4 6 5; 4 6 6; 5 2 6; 5 3 6 ;5 4 6 ; 5 5 6; 5 6 6  ]
% Q = [2 2 3; 2 2 4; 2 3 3; 2 3 4; 2 4 3; 2 4 4; 3 2 4; 3 3 4; 3 4 4; ];
%    Q = [2 2 3; 2 2 4; 2 2 5; 2 3 3; 2 3 4; 2 3 5;  2 4 3; 2 4 4; 2 4 5;  2 5 3; 2 5 4; 2 5 5; 3 2 4; 3 2 5;  3 3 4; 3 3 5;  3 4 4; 3 4 5;  3 5 4; 3 5 5; 4 2 5; 4 3 5; 4 4 5; 4 5 5 ; ];
    Q = [2 2 3; 2 2 4; 2 2 5; 2 2 6; 2 3 3; 2 3 4; 2 3 5; 2 3 6; 2 4 3; 2 4 4; 2 4 5; 2 4 6; 2 5 3; 2 5 4; 2 5 5; 2 5 6; 2 6 3; 2 6 4; 2 6 5; 2 6 6; 3 2 4; 3 2 5; 3 2 6; 3 3 4; 3 3 5; 3 3 6; 3 4 4; 3 4 5; 3 4 6; 3 5 4; 3 5 5; 3 5 6; 3 6 4; 3 6 5; 3 6 6; 4 2 5; 4 2 6; 4 3 5; 4 3 6; 4 4 5; 4 4 6; 4 5 5 ; 4 5 6; 4 6 5; 4 6 6; 5 2 6; 5 3 6 ;5 4 6 ; 5 5 6; 5 6 6  ]

p = Q(:,1);
q = Q(:,2);
r = Q(:,3);
  load('EulerBW4.mat')
F = E;
for i = 1:0
   for j = 1:50
%      tic
     err(i,1) = errordriver(N,p(j),q(j),r(j),1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',4)
     err(i,2) = errordriver(2*N,p(j),q(j),r(j),1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',4)
     E(i,j) = log(err(i,1)/err(i,2))/log(2);
%      toc
     [i j]
   end
    
end
%  save('EulerBW4.mat','E')
%  fit = min(p+q,min(r,p-1))+1+max((r-p).*(q==r),(p==r-1));
%  fit =  p+(r-p).*(q==r);
figure
% plot(fit,'*-')

hold on

boxplot(E)

set(gca,'FontSize',16)
 q = [223; 224; 225; 226; 233; 234; 235; 236; 243; 244; 245; 246; 253; 254; 255; 256; 263; 264; 265; 266; 324; 325; 326; 334; 335; 336; 344; 345; 346; 354; 355; 356; 364; 365; 366; 425; 426; 435; 436; 445; 446; 455 ; 456; 465; 466; 526; 536 ;546 ;  556; 566  ];
B = ['(2,2,3)'; '(2,2,4)'; '(2,2,5)'; '(2,2,6)'; '(2,3,3)'; '(2,3,4)'; '(2,3,5)'; '(2,3,6)'; '(2,4,3)'; '(2,4,4)'; '(2,4,5)'; '(2,4,6)'; '(2,5,3)'; '(2,5,4)'; '(2,5,5)'; '(2,5,6)'; '(2,6,3)'; '(2,6,4)'; '(2,6,5)'; '(2,6,6)'; '(3,2,4)'; '(3,2,5)'; '(3,2,6)'; '(3,3,4)'; '(3,3,5)'; '(3,3,6)'; '(3,4,4)'; '(3,4,5)'; '(3,4,6)'; '(3,5,4)'; '(3,5,5)'; '(3,5,6)'; '(3,6,4)'; '(3,6,5)'; '(3,6,6)'; '(4,2,5)'; '(4,2,6)'; '(4,3,5)'; '(4,3,6)'; '(4,4,5)'; '(4,4,6)'; '(4,5,5)' ; '(4,5,6)'; '(4,6,5)'; '(4,6,6)'; '(5,2,6)'; '(5,3,6)' ;'(5,4,6)' ; '(5,5,6)'; '(5,6,6)'  ]
% B = ['(2,2,3)'; '(2,2,4)';  '(2,3,3)'; '(2,3,4)';  '(2,4,3)'; '(2,4,4)';  '(3,2,4)'; '(3,3,4)';  '(3,4,4)'; ];
%  te = strtrim(cellstr(num2str(q)));
te = cell(size(B,1),1)
% % xlabel(text)
 set(gca,'XTick',1:1:size(B,1),'XTickLabels',te)
 set(gca,'FontSize',20)
ylabel('Order of accuracy of error estimate over 10 runs')



B(:,:)

Xt = 1:size(B,1);
ax = axis;    % Current axis limits
axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
Yl = ax(3:4);

t = text(Xt,Yl(1)*ones(1,length(Xt)),B(:,:));%months(1:2:12,:));
set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
      'Rotation',45,'FontSize',18);
legend('predicted order')
% Remove the default labels
% set(gca,'XTickLabel','')



