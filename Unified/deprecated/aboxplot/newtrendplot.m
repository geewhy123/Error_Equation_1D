close all

x1 = E(:,1);
x2 = E(:,5);
x3 = E(:,9);
x4 = E(:,13);
x5 = E(:,17);

y1 = E(:,2);
y2 = E(:,6);
y3 = E(:,10);
y4 = E(:,14);
y5 = E(:,18);

z1 = E(:,3);
z2 = E(:,7);
z3 = E(:,11);
z4 = E(:,15);
z5 = E(:,19);

w1 = E(:,4);
w2 = E(:,8);
w3 = E(:,12);
w4 = E(:,16);
w5 = E(:,20);
% Concatenate the data sets from each group in a 10000 x 3 matrix
x = cat(2,x1,x2,x3,x4,x5); 
y = cat(2,y1,y2,y3,y4,y5);
z = cat(2,z1,z2,z3,z4,z5);
w = cat(2,w1,w2,w3,w4,w5);
g=figure

subplot(5,2,[1 3 ])
set(gca,'FontSize',30)
set(gca,'YTick',[ 1 2 3 4 5 6],'ylim',[0 6.5])
% Concatenate the each group in a  3 x 10000 x 3 matrix
h = cat(1, reshape(x,[1 size(x)]), reshape(y,[1 size(y)]), reshape(z,[1 size(z)]),reshape(w,[1 size(w)]));
  Q = [2 2 3; 2 2 4; 2 2 5; 2 2 6; 2 3 3; 2 3 4; 2 3 5; 2 3 6; 2 4 3; 2 4 4; 2 4 5; 2 4 6; 2 5 3; 2 5 4; 2 5 5; 2 5 6; 2 6 3; 2 6 4; 2 6 5; 2 6 6; ];
p = Q(:,1);
q = Q(:,2);
r = Q(:,3);

fitt = min(p+q,min(r,p-1))+1+max((r-p).*(q==r),(p==r-1));
colors = [ 0.6000    0.8000    1.0000;  0.4000    0.6000    0.8667;  0.2000    0.4000    0.7333;   0    0.2000    0.6000];

aboxplot(h,'labels',[2,3,4,5,6],'fitt',fitt,'Colormap',   colors)%,'Colormap',[0 0 0;1 0 0;1 0 1;153/255 76/255 0]); % Advanced box plot
% hold on
% plot(fitt,'*')



h = legend('$r=3$','$r=4$','$r=5$','$r=6$','Predicted Order','Orientation','Horizontal'); % Add a legend
set(h,'Interpreter','Latex','FontSize',24)


title('$p=2$','Interpreter','Latex')

subplot(5,2,[2 4 ])

x1 = E(:,21);
x2 = E(:,24);
x3 = E(:,27);
x4 = E(:,30);
x5 = E(:,33);

y1 = E(:,22);
y2 = E(:,25);
y3 = E(:,28);
y4 = E(:,31);
y5 = E(:,34);

z1 = E(:,23);
z2 = E(:,26);
z3 = E(:,29);
z4 = E(:,32);
z5 = E(:,35);


x = cat(2,x1,x2,x3,x4,x5); 
y = cat(2,y1,y2,y3,y4,y5);
z = cat(2,z1,z2,z3,z4,z5);

  Q = [3 2 4; 3 2 5; 3 2 6; 3 3 4; 3 3 5; 3 3 6; 3 4 4; 3 4 5; 3 4 6; 3 5 4; 3 5 5; 3 5 6; 3 6 4; 3 6 5; 3 6 6;] ;

p = Q(:,1);
q = Q(:,2);
r = Q(:,3);

fitt = min(p+q,min(r,p-1))+1+max((r-p).*(q==r),(p==r-1))

set(gca,'FontSize',30)
set(gca,'YTick',[ 1 2 3 4 5 6],'ylim',[0 6.5])
h = cat(1, reshape(x,[1 size(x)]), reshape(y,[1 size(y)]), reshape(z,[1 size(z)]));

aboxplot(h,'labels',[2,3,4,5,6],'fitt',fitt,'Colormap',colors(2:end,:))%,'Colormap',[0 0 0;1 0 0;1 0 1]); % Advanced box plot


title('$p=3$','Interpreter','Latex')

subplot(5,2,[7 9])

x1 = E(:,36);
x2 = E(:,38);
x3 = E(:,40);
x4 = E(:,42);
x5 = E(:,44);

y1 = E(:,37);
y2 = E(:,39);
y3 = E(:,41);
y4 = E(:,43);
y5 = E(:,45);



x = cat(2,x1,x2,x3,x4,x5); 
y = cat(2,y1,y2,y3,y4,y5);

  Q = [4 2 5; 4 2 6; 4 3 5; 4 3 6; 4 4 5; 4 4 6; 4 5 5 ; 4 5 6; 4 6 5; 4 6 6; ];
% 5 2 6; 5 3 6 ;5 4 6 ; 5 5 6; 5 6 6 
p = Q(:,1);
q = Q(:,2);
r = Q(:,3);

fitt = min(p+q,min(r,p-1))+1+max((r-p).*(q==r),(p==r-1));

set(gca,'FontSize',30)
set(gca,'YTick',[ 1 2 3 4 5 6],'ylim',[0 6.5])
h = cat(1, reshape(x,[1 size(x)]), reshape(y,[1 size(y)]));

aboxplot(h,'labels',[2,3,4,5,6],'fitt',fitt,'Colormap',colors(3:end,:))%,'Colormap',[0 0 0;1 0 0]); % Advanced box plot


% ylabel('Order of accuracy of error estimate, n=500','FontSize',25)


title('$p=4$','Interpreter','Latex')
subplot(5,2,[8 10])

x1 = E(:,46);
x2 = E(:,47);
x3 = E(:,48);
x4 = E(:,49);
x5 = E(:,50);


x = cat(2,x1,x2,x3,x4,x5); 


  Q= [5 2 6; 5 3 6 ;5 4 6 ; 5 5 6; 5 6 6 ];
p = Q(:,1);
q = Q(:,2);
r = Q(:,3);

fitt = min(p+q,min(r,p-1))+1+max((r-p).*(q==r),(p==r-1));

set(gca,'FontSize',30)
set(gca,'YTick',[ 1 2 3 4 5 6],'ylim',[0 6.5])
h = cat(1, reshape(x,[1 size(x)]));

aboxplot(h,'labels',[2,3,4,5,6],'fitt',fitt,'Colormap',colors(4:end,:))%,'Colormap',[0 0 0]); % Advanced box plot

annotation(g,'textbox',...
    [0.506208333333333 0.0348869634744324 0.0296875 0.0304449648711944],...
    'Interpreter','latex',...
    'String',{'$$q$$'},...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','AlArabiya',...
    'EdgeColor',[1 1 1]);

ht = text(-6.5,0.85,'Order of accuracy of error estimate, n=500');
set(ht,'Rotation',90,'FontSize',25)

title('$p=5$','Interpreter','Latex')