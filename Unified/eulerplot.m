


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
% B = zeros(50,2);
% for i = 13:24
%     B(i,1) = errordriver(20,p(i),q(i),r(i),1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',4)
%     B(i,2) = errordriver(40,p(i),q(i),r(i),1/3,'D',[1 1],'D',[0.97 1],10,7,'EulerQ','SS',0,'!HC',4)
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