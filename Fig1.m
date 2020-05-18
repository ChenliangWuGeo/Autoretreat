%time steps of model
noTimeStepsFig1 = 5000; %assume each step is one year

%number of grid nodes of model
noNodes = 201; %assmume each step is a grid cell length

%subsample results
spaceJump = 1; %jump every-other grid cell
timeJumpfig1 = 500; % jump 1000 years
tPreAuto = round(tAutobreak/timeJumpfig1+1);
tPostAuto = round(tAutobreak/timeJumpfig1+1)+1;
xJ = 1:spaceJump:noNodes;
tJ = [1:timeJumpfig1:noTimeStepsFig1, noTimeStepsFig1];
[~,a] = size(tJ);
[~,b] = size(xJ);
TJ =nan(a,b);

for i = 1:spaceJump:noNodes
TJ(:,i) = tJ';
end
z = tJ*0+200;

resultSpaceIdx = [1:spaceJump:noNodes];
resultTimeIdx = [1:timeJumpfig1:noTimeStepsFig1, noTimeStepsFig1];

%subsample node location
subSampledX = results.nodeLocation(resultTimeIdx,:);
subSampledX = subSampledX(:,resultSpaceIdx);
autobreakx = results.nodeLocation(tAutobreak,:);

%pull out bed elevation (BE) time step and grid nodes
subSampledBE = results.channelElevation(resultTimeIdx,:);
subSampledBE = subSampledBE(:,resultSpaceIdx);
autobreakBE = results.channelElevation(tAutobreak,:);
        
%pull out desired flow depth time-steps and grid nodes
subSampledflowDepth = results.flowDepth(resultTimeIdx,:);
subSampledflowDepth = subSampledflowDepth(:,resultSpaceIdx);
% autobreakBE = results.channelElevation(resultTimeIdx,tAutobreak);

%pull out desired grain size time-steps and grid nodes
subSampledgrainSize = results.grainSize(resultTimeIdx,:);
subSampledgrainSize = subSampledgrainSize(:,resultSpaceIdx);

%pull out desired shoreline location time-steps and grid nodes
subSampledshorelineLocation2 = results.shorelineLocation2(resultTimeIdx);
subSampledshorelineIdx = results.shorelineIdx(resultTimeIdx);

%pull out desired backwater length time-steps and grid nodes
subSampledbackwaterLength = results.backwaterLength(resultTimeIdx);

%pull out desired basement/foreset break location and elevation
subSampledssbx = results.basementPosition(1,resultTimeIdx);
subSampledssby = results.basementPosition(2,resultTimeIdx);
autobreakssb = results.basementPosition(:,tAutobreak);

%axis and subplot positioning

pos1 = [.1 .66 .35 .26];%0.26 height per figure
pos2 = [.1 .38 .35 .26];
pos3 = [.1 .1 .35 .26];
pos4 = [.55 .66 .35 .26];
pos5 = [.55 .38 .35 .26];
pos6 = [.55 .1 .35 .26];
%plotting

figure(1);hold on
%%
sub(1) = subplot(3,2,1);hold on; 
%plot bed profile, water surface profile and foreset profile
plot(subSampledX(1,:)'./1e3,subSampledBE(1,:)','k','linewidth',1)
plot(subSampledX(1,:)./1e3,subSampledBE(1,:)+subSampledflowDepth(1,:),'b','linewidth',1)
plot([subSampledX(1,:)';subSampledssbx(1)]./1e3,[subSampledBE(1,:)';subSampledssby(1)],'color',[0.55 0.55 0.55],'linewidth',1)

%plot initial sealevel and end sealevel
plot([subSampledX(1,end)./1e3, 150],[subSampledBE(1,end)+subSampledflowDepth(1,end),subSampledBE(1,end)+subSampledflowDepth(1,end)],':b','linewidth',1.5);

%plot shoreline location
plot(subSampledX(1,end)/1e3,subSampledBE(1,end),...
    'ok','linewidth',1,'markerfacecolor','k','markersize',4);
%extrapolate basement profile to 0 km
ssbxExtra = 0;
ssbyExtra = subSampledssby(1) + subSampledssbx(1)*1e-3;

%plot subaqueous basement
plot([ssbxExtra, subSampledssbx]/1e3,[ssbyExtra,subSampledssby],'k','linewidth',1.5);

xlim([-200 50]);
ylim([-40 60]);
ylabel('elevation (m)');
set(gca,'xticklabel',{[]});
set(gca, 'XAxisLocation', 'top')
xticks([0-4*L_auto,0-3.5*L_auto,0-3*L_auto,0-2.5*L_auto,0-2*L_auto,0-1.5*L_auto,0-1*L_auto,0-0.5*L_auto,0,0+0.5* L_auto]);
xticklabels({'-4','-3.5','-3','-2.5','-2','-1.5','-1','-0.5','0','0.5'});
xlabel('{\itx*} ({\itx/L_a})');
    
text(0.1,0.85,'A','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);
text(0.024,0.40,'alluvial profile','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);
text(0.15,0.23,'deltaic foreset','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);
text(0.2,0.03,{'subaqueous','basement'},'Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);
annotation('textarrow',[0.15 0.2],[0.78 0.81],'string','','FontSize',8,'linewidth',1);
annotation('textarrow',[0.32 0.38],[0.73 0.73],'string','','FontSize',8,'linewidth',1);
annotation('textarrow',[0.32 0.38],[0.7 0.7],'string','','FontSize',8,'linewidth',1);
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

ax2 = axes('Position',pos1,'XAxisLocation','bottom','YAxisLocation','right','Color','none');
xlim([-200 50]);
xticks([-200,-150,-100,-50,0]);
xticklabels({});%xticklabels({'-1.5','-1','-0.5','0'});
yticklabels({});
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

%%
sub(2) = subplot(3,2,3);hold on; 
%plot bed profile, water surface profile and foreset profile
plot(subSampledX(1:tPreAuto,:)'./1e3,subSampledBE(1:tPreAuto,:)','k','linewidth',1)
plot(subSampledX(tPreAuto,:)./1e3,subSampledBE(tPreAuto,:)+subSampledflowDepth(tPreAuto,:),'b','linewidth',1)
plot([subSampledX(1:tPreAuto,:)';subSampledssbx(1:tPreAuto)]./1e3,[subSampledBE(1:tPreAuto,:)';subSampledssby(1:tPreAuto)],'color',[0.55 0.55 0.55],'linewidth',1)

%plot initial sealevel and end sealevel
plot([subSampledX(tPreAuto,end)./1e3, 150],[subSampledBE(tPreAuto,end)+subSampledflowDepth(tPreAuto,end),subSampledBE(tPreAuto,end)+subSampledflowDepth(tPreAuto,end)],':b','linewidth',1.5);

%plot autobreak bed profile, water surface profile and foreset profile
plot(autobreakx./1e3,autobreakBE,':k','linewidth',1.5)
plot([autobreakx,autobreakssb(1)]/1e3,[autobreakBE,autobreakssb(2)],':k','linewidth',1)

%plot shoreline location
plot(autobreakx(end)/1e3,autobreakBE(end),...
    'ok','linewidth',1,'markerfacecolor','k','markersize',4);

%extrapolate basement profile to 0 km
ssbxExtra = 0;
ssbyExtra = subSampledssby(1) + subSampledssbx(1)*1e-3;

%plot subaqueous basement
plot([ssbxExtra, subSampledssbx]/1e3,[ssbyExtra,subSampledssby],'k','linewidth',1.5);

xlim([-200 50]);
ylim([-40 60]);
ylabel('elevation (m)');
set(gca,'xticklabel',{[]});
set(gca, 'XAxisLocation', 'top')
xticks([0-4*L_auto,0-3.5*L_auto,0-3*L_auto,0-2.5*L_auto,0-2*L_auto,0-1.5*L_auto,0-1*L_auto,0-0.5*L_auto,0,0+0.5* L_auto]);
text(0.1,0.85,'B','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

ax2 = axes('Position',pos2,'XAxisLocation','bottom','YAxisLocation','right','Color','none');
xlim([-200 50]);
xticks([-200,-150,-100,-50,0]);
xticklabels({});%xticklabels({'-1.5','-1','-0.5','0'});
yticklabels({});
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

%%
sub(3) = subplot(3,2,5);hold on;
p1 = plot([subSampledX(1:tPreAuto,:)';subSampledssbx(1:tPreAuto)]./1e3,[subSampledBE(1:tPreAuto,:)';subSampledssby(1:tPreAuto)],'color',[0.55 0.55 0.55],'linewidth',1);
p5 = plot(subSampledX(end,1:subSampledshorelineIdx(end))./1e3,subSampledBE(end,1:subSampledshorelineIdx(end))+subSampledflowDepth(end,1:subSampledshorelineIdx(end)),'b','linewidth',1);
p3 = plot(subSampledX(tPostAuto:end,:)'./1e3,subSampledBE(tPostAuto:end,:)','k','linewidth',1);
plot([subSampledX(tPostAuto:end,:)';subSampledssbx(tPostAuto:end)]./1e3,[subSampledBE(tPostAuto:end,:)';subSampledssby(tPostAuto:end)],'k','linewidth',1)
plot(autobreakx./1e3,autobreakBE,':k','linewidth',1.5)
p2 = plot([autobreakx,autobreakssb(1)]/1e3,[autobreakBE,autobreakssb(2)],':k','linewidth',1);

%plot shoreline location
p4 = plot(results.shorelineLocation2(5000)/1e3,subSampledBE(end,results.shorelineIdx(5000)),...
    'ok','linewidth',1,'markerfacecolor','k','markersize',4);
plot([ssbxExtra, subSampledssbx]/1e3,[ssbyExtra,subSampledssby],'k','linewidth',1.5);

% plot([-500 500],[0 0],':b','linewidth',1.5)
p6 = plot([subSampledX(end,subSampledshorelineIdx(end))./1e3, 150],[subSampledBE(end,subSampledshorelineIdx(end))+subSampledflowDepth(end,subSampledshorelineIdx(end)),subSampledBE(end,subSampledshorelineIdx(end))+subSampledflowDepth(end,subSampledshorelineIdx(end))],':b','linewidth',1.5);
xlim([-200 50]);
ylim([-40 60]);

% set(gca,'xticklabel',{[]});
set(gca, 'XAxisLocation', 'top')
xticks([0-4*L_auto,0-3.5*L_auto,0-3*L_auto,0-2.5*L_auto,0-2*L_auto,0-1.5*L_auto,0-1*L_auto,0-0.5*L_auto,0,0+0.5* L_auto]);
xticklabels({});
ylabel('elevation (m)');
text(0.1,0.85,'C','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);

L = legend([p4(1),p1(1) p2(1), p3(1),p5(1),p6(1)],...
    {'shoreline','pre-autobreak','autobreak','post-autobreak',...
    'flow surface profile','base level'},'location','southwest');
legend('boxoff') 
legendPosition = [0.13, 0.085, .22, .16];
set(L, 'Position', legendPosition)
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

ax2 = axes('Position',pos3,'XAxisLocation','bottom','YAxisLocation','right','Color','none');
xlim([-200 50]);
xticks([-200,-150,-100,-50,0]);
xlabel('distance {\itx} (km)');
yticklabels({});
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);


subSampledX2 = nan(11,100);
subSampledH2 = nan(11,100);
subSampledD2 = nan(11,100);
TJ2 = TJ(:,1:100);

for i = 1:11
    idx = subSampledshorelineIdx(i);
    subSampledX2(i,:) = subSampledX(i,idx-99:idx);
    subSampledH2(i,:) = subSampledflowDepth(i,idx-99:idx);
    subSampledD2(i,:) = subSampledgrainSize(i,idx-99:idx);
end

%%
sub(4) = subplot(3,2,2);hold on; 
%plot bed profile, water surface profile and foreset profile
shorelineElevation = nan(1,5000);
for i = 1:5000
shorelineElevation(i) = results.channelElevation(i,results.shorelineIdx(i));
end
plot(results.shorelineLocation2(1:5000)./1e3,shorelineElevation,'k','linewidth',1)

xlim([-200 50]);
ylim([-10 40]);
ylabel('elevation (m)');
set(gca,'xticklabel',{[]});
set(gca, 'XAxisLocation', 'top')
xticks([0-4*L_auto,0-3.5*L_auto,0-3*L_auto,0-2.5*L_auto,0-2*L_auto,0-1.5*L_auto,0-1*L_auto,0-0.5*L_auto,0,0+0.5* L_auto]);
xticklabels({'-4','-3.5','-3','-2.5','-2','-1.5','-1','-0.5','0','0.5'});
xlabel('{\itx*} ({\itx/L_a})');
    
text(0.1,0.85,'D','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);
text(0.05,0.4,{'sediment-','starved autoretreat'},'Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);
text(0.55,0.27,'retreat','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);
text(0.45,0.11,{'advance'},'Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);

annotation('textarrow',[0.83 0.70],[0.80 0.90],'string','','FontSize',8,'linewidth',1);
annotation('textarrow',[0.85 0.84],[0.72 0.76],'string','','FontSize',8,'linewidth',1);
annotation('textarrow',[0.84 0.85],[0.67 0.71],'string','','FontSize',8,'linewidth',1);

set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

ax2 = axes('Position',pos4,'XAxisLocation','bottom','YAxisLocation','right','Color','none');
xlim([-200 50]);
xticks([-200,-150,-100,-50,0]);
xticklabels({});%xticklabels({'-1.5','-1','-0.5','0'});
yticklabels({});
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

%%
figure(1);hold on
% grain size in stratigraphy
sub(5) = subplot(3,2,4);hold on
surf(subSampledX./1e3,subSampledBE,subSampledgrainSize/subSampledgrainSize(1),'EdgeColor','interp','FaceColor','interp');

% initial and final profiles
plot(subSampledX(1,:)'./1e3,subSampledBE(1,:)','k','linewidth',1)
plot([subSampledX(1,:)';subSampledssbx(1)]./1e3,[subSampledBE(1,:)';subSampledssby(1)],'color',[0.55 0.55 0.55],'linewidth',1)
plot(subSampledX(end,:)'./1e3,subSampledBE(end,:)','k','linewidth',1);
plot([subSampledX(end,:)';subSampledssbx(end)]./1e3,[subSampledBE(end,:)';subSampledssby(end)],'k','linewidth',1)
plot([ssbxExtra, subSampledssbx]/1e3,[ssbyExtra,subSampledssby],'k','linewidth',1.5);

set(gca,'ColorScale','log')
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);
view(2)
shading interp
xlim([-200 50]);
ylim([-40 60]);
xticks([-200,-150,-100,-50,0]);
xticklabels({});
ylabel('elevation (m)');
text(0.1,0.85,'E','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);
[~,temp] = size(subSampledX(end,subSampledshorelineIdx(end):end));
temp = ones(1,temp) + 200;
plot3(subSampledX(end,subSampledshorelineIdx(end):end)./1e3, subSampledBE(end,subSampledshorelineIdx(end):end),temp,'r','LineWidth',1)

ax2 = axes('Position',pos5,'XAxisLocation','top','YAxisLocation','right','Color','none');
xlim([-200 50]);
ylim([-10 50]);
xticks([-1.5*L_auto,-1*L_auto,-0.5*L_auto,0]);
xticklabels({});%xticklabels({'-1.5','-1','-0.5','0'});
yticklabels({});% yticklabels({'0','0.5','1','1.5','2','2.5'});
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

ax3 = axes('Position',[.63 .43 .14 .045],'XAxisLocation','bottom','YAxisLocation','right','Color','none');
plot(subSampledX(end,subSampledshorelineIdx(end):end)./1e3, subSampledgrainSize(end,subSampledshorelineIdx(end):end)/subSampledgrainSize(1),'r','LineWidth',1) 
xlim([-150 10]);
ylim([.2 .3]);
yticks([0.2 0.3]);
xticks([-100 0]);
xlabel('{\itx} (km)');
ylabel('\it{D*}')
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

h5 = colorbar;
h5.Label.String = '\it{D*}';
h5.Label.Rotation = 90;
h5.TicksMode = 'manual';
h5.Ticks = [0.3, 0.6, 1];
h5.TickLabels = {'0.3','0.6','1'} ;
set(h5, 'YAxisLocation','left')
set(h5,'position',[0.85 0.51 .02 0.12]);
[0.85 0.23 .02 0.12];
%%
figure(1);hold on
sub(6) = subplot(3,2,6);hold on
surf(subSampledX2./1e3,tJ,subSampledH2/subSampledflowDepth(1),'EdgeColor','interp','FaceColor','interp')
[c, h] = contour3(subSampledX2./1e3,TJ2,subSampledH2/subSampledflowDepth(1),'k');
view(2)
clabel(c,h,[1.2, 1.6, 2]);
xlim([-200 50]);
ylabel('time t (yr)');
xlabel('distance {\itx} (km)');
shading interp
ax1 = gca;
text(0.1,0.85,'F','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);
plot3((subSampledshorelineLocation2-subSampledbackwaterLength)/1e3,tJ,z,'w','linewidth',1);
plot3(subSampledshorelineLocation2(tPreAuto:end)./1e3, tJ(tPreAuto:end),z(tPreAuto:end),'r','LineWidth',1.5)
xlim([-200 50]);
ylim([0 5000]);
xticks([-200,-150,-100,-50,0]);
yticks([0 1000 2000 3000 4000 5000]);
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

ax2 = axes('Position',pos6,'XAxisLocation','top','YAxisLocation','right','Color','none');
xlim([-200 50]);
ylim([0 5000]);
xticks([-1.5*L_auto,-1*L_auto,-0.5*L_auto,0]);
xticklabels({});%xticklabels({'-1.5','-1','-0.5','0'});
yticks([0, 0.5*T_auto,1*T_auto,1.5*T_auto,2*T_auto]);
yticklabels({'0','0.5','1','1.5','2','2.5'});
ylabel('{\itt*} ({\itt/T_a})');
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

h6 = colorbar(ax1);
h6.Label.String = '\it{H*}';
h6.Label.Rotation = 90;
h6.TicksMode = 'manual';
h6.Ticks = [1, 1.5, 2];
h6.TickLabels = {'1','1.5','2'} ;
set(h6, 'YAxisLocation','left')
set(h6,'position',[0.85 0.23 .02 0.12]);

% linkaxes([ax1,ax2])
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);
view(2)

set(sub(1),'position',pos1)
set(sub(2),'position',pos2)
set(sub(3),'position',pos3)
set(sub(4),'position',pos4)
set(sub(5),'position',pos5)
set(sub(6),'position',pos6)

set(gcf,'Position',[10 10 500 700])
