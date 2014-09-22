


%osc 223 226

%plot 1

close all
load('euler244.mat')

% figure

h=figure
set(gca,'FontSize',18)
F = 18;
        subplot(3,2,1)
        plot(x,eu(:,1),'*',x,exacterru(:,1),'o')
        xlabel('$\epsilon_\rho$','Interpreter','Latex','FontSize',F);
        title('Exact and Computed Error','Interpreter','Latex','FontSize',F,'FontWeight','bold')
                l=legend('Error estimate','Exact error')
                set(l,'FontSize',F)
                
        subplot(3,2,3)
        plot(x,eu(:,2),'*',x,exacterru(:,2),'o')
        xlabel('$\epsilon_{\rho u}$','Interpreter','Latex');
        subplot(3,2,5)
        plot(x,eu(:,3),'*',x,exacterru(:,3),'o')
        xlabel('$\epsilon_{\rho E}$','Interpreter','Latex');
       


        subplot(3,2,2)
        plot(x,exacterru(:,1)-eu(:,1),'x')
        xlabel('$\epsilon_\rho$','Interpreter','Latex');

        subplot(3,2,4)
        plot(x,exacterru(:,2)-eu(:,2),'x')
        xlabel('(computed)$\epsilon_u$','Interpreter','Latex');
        subplot(3,2,6)
        plot(x,exacterru(:,3)-eu(:,3),'x')
        xlabel('(computed)$\epsilon_P$','Interpreter','Latex');
        


error('1')

%trend


% close all
clear all
% load('Uunif500.mat')
Q = [2 2 3; 2 2 4; 2 2 5; 2 3 3; 2 3 4; 2 3 5; 2 4 3; 2 4 4; 2 4 5; 2 5 3; 2 5 4; 2 5 5; 3 2 4; 3 2 5; 3 3 4; 3 3 5;3 4 4; 3 4 5; 3 5 4; 3 5 5; 4 2 5; 4 3 5; 4 4 5; 4 4 5];

p = Q(:,1);
q = Q(:,2);
r = Q(:,3);
tic
load('Uunif25.mat')
% for i = 1:25
%    for j = 1:24
%      
%      err(i,1) = errordriver(20,p(j),q(j),r(j),1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',4)
%      err(i,2) = errordriver(40,p(j),q(j),r(j),1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',4)
%      E(i,j) = log(err(i,1)/err(i,2))/log(2);
%      
%      
%    end
%     
% end
toc
% fit = min(p+q,min(r,p-1))+1+max((r-p).*(q==r),(p==r-1));

% plot(fit,'*-')
figure
hold on

boxplot(E)

set(gca,'FontSize',16)
 q = [223; 224; 225; 226; 233; 234; 235; 236; 243; 244; 245; 246; 253; 254; 255; 256; 263; 264; 265; 266; 324; 325; 326; 334; 335; 336; 344; 345; 346; 354; 355; 356; 364; 365; 366; 425; 426; 435; 436; 445; 446; 455 ; 456; 465; 466; 526; 536 ;546 ;  556; 566  ];
B = ['(2,2,3)'; '(2,2,4)'; '(2,2,5)';  '(2,3,3)'; '(2,3,4)'; '(2,3,5)';  '(2,4,3)'; '(2,4,4)'; '(2,4,5)';  '(2,5,3)'; '(2,5,4)'; '(2,5,5)'; '(3,2,4)'; '(3,2,5)';  '(3,3,4)'; '(3,3,5)'; '(3,4,4)'; '(3,4,5)'; '(3,5,4)'; '(3,5,5)'; '(4,2,5)';  '(4,3,5)'; '(4,4,5)';  '(4,5,5)' ;  ]

%  te = strtrim(cellstr(num2str(q)));
te = cell(50,1)
% % xlabel(text)
 set(gca,'XTick',1:1:50,'XTickLabels',te)
 set(gca,'FontSize',16)
ylabel('Order of accuracy of error estimate over 25 runs')



B(:,:)

Xt = 1:24;
ax = axis;    % Current axis limits
axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
Yl = ax(3:4);

t = text(Xt,Yl(1)*ones(1,length(Xt)),B(:,:));%months(1:2:12,:));
set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
      'Rotation',45,'FontSize',14);
legend('predicted order')
% Remove the default labels
% set(gca,'XTickLabel','')


error('1')

% A(1,1) = errordriver(10,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'HC',1)
% A(1,2) = errordriver(20,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'HC',1)
% A(1,3) = errordriver(40,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'HC',1)
% A(2,1) = errordriver(10,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'HC',2)
% A(2,2) = errordriver(20,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'HC',2)
% A(2,3) = errordriver(40,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'HC',2)
% A(3,1) = errordriver(10,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'HC',3)
% A(3,2) = errordriver(20,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'HC',3)
% A(3,3) = errordriver(40,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'HC',3)
% A(4,1) = errordriver(10,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',4)
% A(4,2) = errordriver(20,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',4)
% A(4,3) = errordriver(40,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',4)
% A(5,1) = errordriver(10,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',2)
% A(5,2) = errordriver(20,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',2)
% A(5,3) = errordriver(40,2,4,4,1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',2)
close all
%  load('euler.mat')
x = [1/10 1/20 1/40];
xmin = 1e-2;
xmax = 2e-1;
y = linspace(xmin,xmax,100);
c = 1e-5;
b4 = c*y.^4;
b2 = c*y.^2;
figure
set(gca,'FontSize',14)
loglog(x,A,'*-',y,b2,y,b4,'LineWidth',1.5)
xlim([xmin xmax])
xlabel('r')
ylabel('max ||e-e_h||_2')
legend('method 1','method 2','method 3','method 4','method 2 correct','2nd order','4th order')


Q = [2 2 3; 2 2 4; 2 2 5; 2 3 3; 2 3 4; 2 3 5; 2 4 3; 2 4 4; 2 4 5; 2 5 3; 2 5 4; 2 5 5; 3 2 4; 3 2 5; 3 3 4; 3 3 5;3 4 4; 3 4 5; 3 5 4; 3 5 5; 4 2 5; 4 3 5; 4 4 5; 4 4 5];
p = Q(:,1);
q = Q(:,2);
r = Q(:,3);
% B = zeros(24,2);
% for i = 1:24
%     B(i,1) = errordriver(20,p(i),q(i),r(i),0,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',4)
%     B(i,2) = errordriver(40,p(i),q(i),r(i),0,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',4)
% end
C = [223;224; 225;233; 234; 235; 243; 244; 245; 253; 254; 255;324;325;334; 335;  344; 345;  354; 355; 425; 435; 445; 455 ;  ];
z = 1:24;
figure
plot(z,log(B(1:24,1)./B(1:24,2))./log(2),'*-')

 te = strtrim(cellstr(num2str(C)));
% te = cell(50,1)
% % xlabel(text)
 set(gca,'XTick',1:1:24,'XTickLabels',te)

 ylabel('error order of accuracy 1 run 20/40 rand')
 
 
 
 
 
 
 
 
 