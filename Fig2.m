%time steps of model
noTimeStepsFig1 = 5000; %assume each step is one year

%number of grid nodes of model
noNodes = 201; %assmume each step is a grid cell length

%subsample results
spaceJump = 1; %jump every-other grid cell
timeJumpfig1 = 100; % jump 1000 years
tPreAuto = round(tAutobreak/timeJumpfig1+1);
tPostAuto = round(tAutobreak/timeJumpfig1+1)+1;
xJ = 1:spaceJump:noNodes;
tJ = [1:timeJumpfig1:noTimeStepsFig1, noTimeStepsFig1];
% for i = 1:spaceJump:noNodes
% TJ(:,i) = tJ';
% end
% z = tJ*0+200;

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

%pull out desired backwater length time-steps and grid nodes
subSampledbackwaterLength = results.backwaterLength(resultTimeIdx);

%pull out desired basement/foreset break location and elevation
subSampledssbx = results.basementPosition(1,resultTimeIdx);
subSampledssby = results.basementPosition(2,resultTimeIdx);
autobreakssb = results.basementPosition(:,tAutobreak);

%pull out desired basement/foreset break location and elevation
subSampledsedPartition = results.sedPartition(resultTimeIdx);

%pull out desired max aggradation location 
subSampledmaxAggradationX = results.maxAggradationX(resultTimeIdx);

%pull out desired shoreline migration rate
subSampledshorelineMigRate2 = shorelineMigRate2(resultTimeIdx);

figure(2);hold on

%axis and subplot positioning
pos1 = [.15 .7 .7 .17];
pos2 = [.15 .5 .7 .17];
pos3 = [.15 .3 .7 .17];
pos4 = [.15 .1 .7 .17];


sub1 = subplot(4,1,1);hold on;
plot(tJ,subSampledshorelineMigRate2,'ok','linewidth',1,'markerSize',4);
plot([tAutobreak,tAutobreak],[-60 20],':','color',[0.5 0.5 0.5],'lineWidth',1.5)
plot([tRetreat,tRetreat],[-60 20],':','color',[0.5 0.5 0.5],'lineWidth',1.5)
plot([0 5000],[0 0],':k','lineWidth',1)
ylabel({'migration','rate (m/yr)'});
xlim([0 5000]);
xticks([0 1000 2000 3000 4000 5000]);
xticklabels([]);
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
text(0.9,0.75,'A','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);
text(0.02,-0.2,'advance','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);
text(0.25,-0.2,'retreat','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);
text(0.4,-0.2,'sediment-starved autoretreat','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);

sub2 = subplot(4,1,2);hold on;
p2 = plot(tJ,subSampledmaxAggradationX/1e3,'k','linewidth',1);
p1 = plot(tJ,subSampledshorelineLocation2/1e3,':k','linewidth',1.5); 
plot([tAutobreak,tAutobreak],[-200 100],':','color',[0.5 0.5 0.5],'lineWidth',1.5)
plot([tRetreat,tRetreat],[-200 100],':','color',[0.5 0.5 0.5],'lineWidth',1.5)
ylabel('distance (km)');
xlim([0 5000]);
ylim([-200 50]);
xticklabels({});
xticks([0 1000 2000 3000 4000 5000]);

ax2 = axes('Position',pos2,'XAxisLocation','top','YAxisLocation','right','Color','none');
yticks([0-4*L_auto,0-3.5*L_auto,0-3*L_auto,0-2.5*L_auto,0-2*L_auto,0-1.5*L_auto,0-1*L_auto,0-0.5*L_auto,0,0+0.5* L_auto]);
yticklabels({'-4','-3.5','-3','-2.5','-2','-1.5','-1','-0.5','0','0.5'});
ylabel('\it{x*}');
ylim([-200 50]);
ax2.YAxis.Color = 'k';
xlim([0 5000]);
xticklabels({});
xticks([0, 0.5*T_auto, T_auto, 1.5*T_auto, 2*T_auto,2.5*T_auto]);
h1 = legend([p1 p2],{'shoreline location','depositional front'},'location','southwest');
text(0.9,0.75,'B','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);

sub3 = subplot(4,1,3);hold on;
plot(tJ,subSampledbackwaterLength./1e3,'k','linewidth',1);
plot([tAutobreak,tAutobreak],[0 200],':','color',[0.5 0.5 0.5],'lineWidth',1.5)
plot([tRetreat,tRetreat],[0 200],':','color',[0.5 0.5 0.5],'lineWidth',1.5)
ylabel('{\itL_b} (km)');
ylim([0 210]);
xlim([0 5000]);
xticklabels({});
xticks([0 1000 2000 3000 4000 5000]);
ax2 = axes('Position',pos3,'XAxisLocation','top','YAxisLocation','right','Color','none');
yticks([0,0.5* L_auto,1*L_auto,1.5*L_auto,2*L_auto,2.5*L_auto,3*L_auto,3.5*L_auto,4*L_auto]);
yticklabels({'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5'});
ylabel('\it{x*}');
ylim([0 210]);
ax2.YAxis.Color = 'k';
xlim([0 5000]);
xticklabels({});
xticks([0, 0.5*T_auto, T_auto, 1.5*T_auto, 2*T_auto,2.5*T_auto]);
text(0.9,0.75,'C','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);

sub4 = subplot(4,1,4);hold on;
plot(tJ,subSampledsedPartition*100,'k','linewidth',1);
plot([tAutobreak,tAutobreak],[0 100],':','color',[0.5 0.5 0.5],'lineWidth',1.5)
plot([tRetreat,tRetreat],[0 100],':','color',[0.5 0.5 0.5],'lineWidth',1.5)
ylabel({'sediment','partition (%)'});
xlim([0 5000]);
xticks([0 1000 2000 3000 4000 5000]);
xlabel('time {\itt} (year)');
ax2 = axes('Position',pos4,'XAxisLocation','top','YAxisLocation','right','Color','none');
yticks([]);
yticklabels({});
ax2.YAxis.Color = 'k';
xlim([0 5000]);
xticklabels({});
xticks([0, 0.5*T_auto, T_auto, 1.5*T_auto, 2*T_auto,2.5*T_auto]);
text(0.9,0.75,'D','Units', 'Normalized', 'VerticalAlignment', 'bottom','HorizontalAlignment','left','fontSize',10);

sub1 = subplot(4,1,1);hold on;
ax2 = axes('Position',pos1,'XAxisLocation','top','YAxisLocation','right','Color','none');
xticks([0, 0.5*T_auto, T_auto, 1.5*T_auto, 2*T_auto,2.5*T_auto]);
xticklabels({'0','0.5','1','1.5','2','2.5'});
xlabel('{\itt*} ({\itt/T_a})');
xlim([0 5000]); %if not specify this, it will not plot properly
yticks([]);
yticklabels({});
ax2.XAxis.Color = 'k';

set(sub1,'position',pos1)
set(sub2,'position',pos2)
set(sub3,'position',pos3)
set(sub4,'position',pos4)