load('MC_result1.mat');
load('MC_slopeInput.mat');
load('MC_baseSlopeInput.mat');

fig4Results = MC_results;
%excluding cases where autobreak is not reached.
fig4Results(fig4Results(:,4)<=5, 4) = NaN;

figure(4); hold on
subplot(1,2,1);hold on; box on

p = scatter(fig4Results(:,2)./fig4Results(:,5),fig4Results(:,4)...
    ./fig4Results(:,6),8,baseSlopeInput,'filled');
xlim([1.4, max(fig4Results(:,2)./fig4Results(:,5))]);
colormap(parula);
cb = colorbar;
caxis([0.00001, 0.001]);
set(gca,'ColorScale','log')
set(cb,'position',[.2, .44, .03, .3]);
title(cb, {'basin slope'})
cb.Ticks = [0.00001, 0.0001, 0.001];
cb.TickLabels = {'10^{-5}','10^{-4}','10^{-3}'} ;
set(cb, 'YAxisLocation','left')
xlabel('\itL_b/L_a');
ylabel('\itt_{ssa}/T_a');
text(0.1,0.9, 'A', 'FontSize', 8, 'Color', 'k','unit','normalized','VerticalAlignment','bottom','HorizontalAlignment','left');
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

figure(4);hold on
subplot(1,2,2);hold on; box on
Hnormal = fig4Results(:,7);
Hmouth = fig4Results(:,8);
symb = ['*','o','s','^'];

p = plot(Hnormal,Hmouth,'o');
set(p,'color',[.7 .7 .7],'markerfacecolor',[.7 .7 .7],'markersize',3);

depthFit = fit(Hnormal,Hmouth,'poly1');
f(1) = plot(depthFit,'k');
set(f(1),'LineWidth',1);
f(2) = plot([0,20],[0 40],'-.k');

realExampleLabel = {'Trinity','Mississippi','Fly','Yangtze'...
    'Tombigbee','Brazos'};
realExampleX = [2.32, 10.31, 7.63, 12, 6.57, 2.61];
realExampleXmax = [3.93, 10.99, 9.79, 12, 6.90, 2.61];
realExampleXmin = [1.07, 9.10, 4.68, 12, 6.40, 2.61];
realExampleY = [5.67, 21.2, 14.37, 23, 12.54, 5.95];
realExampleYmax = [6.47, 22.93, 18.97, 23, 13.34, 8.75];
realExampleYmin = [4.92, 20.92, 10.43, 23, 11.73, 4.42];

realExampleXLabel = [2, 10.31, 7.63, 12, 6.57, 3.5];
realExampleYLabel = [5, 21.2, 14.37, 23, 12.54, 8.4];

yneg = realExampleY - realExampleYmin;
ypos = realExampleYmax - realExampleY;
xneg = realExampleX - realExampleXmin;
xpos = realExampleXmax - realExampleX;
e = errorbar(realExampleX,realExampleY,yneg,ypos,xneg,xpos,'o',...
    'MarkerSize',6,'MarkerEdgeColor','white','MarkerFaceColor','red','CapSize',0);
e.Color = 'red';
text(realExampleXLabel,realExampleYLabel,realExampleLabel,'VerticalAlignment','top','HorizontalAlignment','left','fontsize',8)

R = corrcoef(Hnormal,Hmouth);
R2 = R(1,2).^2;
caption = sprintf('y = %.3f * x + %.3f\n r^2 = %.3f', depthFit.p1, depthFit.p2,R2);
text(0.5,0.03, caption, 'FontSize', 8, 'Color', 'k','unit','normalized','VerticalAlignment','bottom','HorizontalAlignment','left');

xlabel('normal flow depth {\itHn} (m)');
ylabel('mouth flow depth {\itHm} (m)');
text(0.15,0.9, 'B', 'FontSize', 8, 'Color', 'k','unit','normalized','VerticalAlignment','bottom','HorizontalAlignment','left');

L1 = legend([f(2),f(1),p,e],{'{\itHm/Hn} = 2','model fit','modeled','natural'},'location','best');
legend('boxoff') 
legendPosition = [0.52, 0.675, .3, .1];
set(L1, 'Position', legendPosition)

set(gcf,'position',[100 100 600 250]);
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);


