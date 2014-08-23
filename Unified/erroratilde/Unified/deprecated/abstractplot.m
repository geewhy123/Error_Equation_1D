% ploterra.m has 20 plots



% % % 
% % % %osc 223 226
% % % 
% % % errordriver(40,2,2,3,1/3,100,.3,7,'Advection')
% % % xlabel('$$x$$','FontSize',38,'Interpreter','Latex')
% % % ylabel('$$\mathbf{\mathcal{I}^h\epsilon - \epsilon_h}$$','FontSize',38,'Interpreter','Latex')
% % % % set(gca,'FontSize',16')
% % % 
% % % error('1')
% % % figure
% % % errordriver(40,2,2,6,1/3,100,.3,7,'Advection')
% % % xlabel('$$x$$','FontSize',38,'Interpreter','Latex')
% % % ylabel('$$\mathbf{\mathcal{I}^h\epsilon - \epsilon_h}$$','FontSize',38,'Interpreter','Latex')
% % % % set(gca,'FontSize',16')
% % % error('1')
% % % %trend
% % % 
% % % 
% % % close all
% % % clear all
% % % load('Uunif500.mat')
% % %   Q = [2 2 3; 2 2 4; 2 2 5; 2 2 6; 2 3 3; 2 3 4; 2 3 5; 2 3 6; 2 4 3; 2 4 4; 2 4 5; 2 4 6; 2 5 3; 2 5 4; 2 5 5; 2 5 6; 2 6 3; 2 6 4; 2 6 5; 2 6 6; 3 2 4; 3 2 5; 3 2 6; 3 3 4; 3 3 5; 3 3 6; 3 4 4; 3 4 5; 3 4 6; 3 5 4; 3 5 5; 3 5 6; 3 6 4; 3 6 5; 3 6 6; 4 2 5; 4 2 6; 4 3 5; 4 3 6; 4 4 5; 4 4 6; 4 5 5 ; 4 5 6; 4 6 5; 4 6 6; 5 2 6; 5 3 6 ;5 4 6 ; 5 5 6; 5 6 6  ]
% % % 
% % % p = Q(:,1);
% % % q = Q(:,2);
% % % r = Q(:,3);
% % % 
% % % fit = min(p+q,min(r,p-1))+1+max((r-p).*(q==r),(p==r-1));
% % % 
% % % plot(fit,'*-')
% % % 
% % % hold on
% % % 
% % % boxplot(E)
% % % 
% % % set(gca,'FontSize',16)
% % %  q = [223; 224; 225; 226; 233; 234; 235; 236; 243; 244; 245; 246; 253; 254; 255; 256; 263; 264; 265; 266; 324; 325; 326; 334; 335; 336; 344; 345; 346; 354; 355; 356; 364; 365; 366; 425; 426; 435; 436; 445; 446; 455 ; 456; 465; 466; 526; 536 ;546 ;  556; 566  ];
% % % B = ['(2,2,3)'; '(2,2,4)'; '(2,2,5)'; '(2,2,6)'; '(2,3,3)'; '(2,3,4)'; '(2,3,5)'; '(2,3,6)'; '(2,4,3)'; '(2,4,4)'; '(2,4,5)'; '(2,4,6)'; '(2,5,3)'; '(2,5,4)'; '(2,5,5)'; '(2,5,6)'; '(2,6,3)'; '(2,6,4)'; '(2,6,5)'; '(2,6,6)'; '(3,2,4)'; '(3,2,5)'; '(3,2,6)'; '(3,3,4)'; '(3,3,5)'; '(3,3,6)'; '(3,4,4)'; '(3,4,5)'; '(3,4,6)'; '(3,5,4)'; '(3,5,5)'; '(3,5,6)'; '(3,6,4)'; '(3,6,5)'; '(3,6,6)'; '(4,2,5)'; '(4,2,6)'; '(4,3,5)'; '(4,3,6)'; '(4,4,5)'; '(4,4,6)'; '(4,5,5)' ; '(4,5,6)'; '(4,6,5)'; '(4,6,6)'; '(5,2,6)'; '(5,3,6)' ;'(5,4,6)' ; '(5,5,6)'; '(5,6,6)'  ]
% % % 
% % % %  te = strtrim(cellstr(num2str(q)));
% % % te = cell(50,1)
% % % % % xlabel(text)
% % %  set(gca,'XTick',1:1:50,'XTickLabels',te)
% % %  set(gca,'FontSize',16)
% % % ylabel('Order of accuracy of error estimate over 500 runs')
% % % 
% % % 
% % % 
% % % B(:,:)
% % % 
% % % Xt = 1:50;
% % % ax = axis;    % Current axis limits
% % % axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
% % % Yl = ax(3:4);
% % % 
% % % t = text(Xt,Yl(1)*ones(1,length(Xt)),B(:,:));%months(1:2:12,:));
% % % set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
% % %       'Rotation',45,'FontSize',14);
% % % legend('predicted order')
% % % % Remove the default labels
% % % % set(gca,'XTickLabel','')
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % %histogram
% % % 
% % % w=figure
% % % 
% % % 
% % % for j = 1:50
% % %     X(j) = chi2gof(E(:,j));
% % % %    if(chi2gof(E(:,j))==1)
% % % %       error('1') 
% % % %    end
% % % 
% % % 
% % % subplot(10,5,j)
% % % % hist(E(:,j))
% % % if(X(j)==1)
% % % set(get(gca,'child'),'FaceColor','r');
% % % end
% % % 
% % %  xlabel(B(j,:),'FontSize',12);%num2str(Q(j,:)))
% % % % if (k==1)
% % % %    legend('exact error','computed error'); 
% % % % end
% % % hold on
% % % 
% % % l=histfit(E(:,j))
% % % 
% % % % t(j)=ttest(E(:,j),fit(j));
% % % x = [fit(j),fit(j)];
% % % y = [0,100];
% % % l2=plot(x,y,'*b-')
% % % 
% % % 
% % % e = E(:,j);
% % % e = (e-mean(e))/sqrt(var(e));
% % % [h(j),p(j)]=kstest(e);
% % % 
% % % if(j==26)
% % %    ylabel('Frequency','FontSize',18) 
% % % end
% % % 
% % % 
% % % end
% % % L = [l(2) l2];
% % % % annotation(w,'textbox',...
% % % %     [0.0504791666666667 0.532837474674436 0.0322916666666667 0.0351288056206089],...
% % % %     'Interpreter','latex',...
% % % %     'String',{'Frequency'},...
% % % %     'FontWeight','bold',...
% % % %     'FontSize',18,...
% % % %     'FontName','AlArabiya',...
% % % %     'EdgeColor',[1 1 1]);
% % % 
% % % annotation(w,'textbox',...
% % %     [0.506208333333333 0.0548869634744324 0.0296875 0.0304449648711944],...
% % %     'Interpreter','latex',...
% % %     'String',{'$$||\mathbf{\epsilon}-\tilde{\bf{\epsilon}}||$$ order'},...
% % %     'FontWeight','bold',...
% % %     'FontSize',18,...
% % %     'FontName','AlArabiya',...
% % %     'EdgeColor',[1 1 1]);
% % % legend(L,'normal distribution','theoretical order')


figure

t2 = [0.149 0.367 1.289 5.39];
t4 = [0.155 0.418 1.391 5.334 21.421];
e2 = [1.85e-3 4.44e-4 1.23e-4 2.93e-5];
e4 = [0.00128 4.12e-5 2.65e-6 1.75e-7 1.06e-8];

loglog(t2,e2,'--^',t4,e4,'--s','LineWidth',3,'MarkerSize',10)
xlabel('average process time, t (s)','FontSize',38,'Interpreter','Latex')
ylabel('error, $${||\epsilon||_2}$$','FontSize',38,'Interpreter','Latex')
legend('2^{nd} order scheme','4^{th} order scheme')
set(gca,'FontSize',40)
