noTimeStepsPlot = 8000;

BEI = results.aggradationRate;
channelElevation = results.channelElevation;
detrendx =  results.nodeLocation;
detrendAdjust = nan(noTimeStepsPlot,201);

for j=1:noTimeStepsPlot
    for i = 2:201
        detrendAdjust(j,i) = channelElevation(j,1) - Sfi *(detrendx(j,i)-detrendx(j,1));
    end
end
detrendAdjust(:,1) = channelElevation(:,1);
detrendCE = channelElevation - detrendAdjust;

[length,~ ] = size(BEI);
xStep = results.nodeLocation;

rMax = nan(1,length);
rMaxCE = nan(1,length);
drMax = nan(1,length);
ave_r = nan(1,length);
d_ave_r = nan(1,length);
RMS = nan(1,length);
RMSce = nan(1,length);
lagMax = nan(1,length);

cmap = parula(length);
for i=1:1:length-1
    rMax(i) = max(xcorr(BEI(i,:),BEI(i+1,:)));
    ave_r(i) = sum((xcorr(BEI(i,:),BEI(i+1,:))) * (xStep(i,end)-xStep(i,1))/100);
    RMS(i) = sqrt(sum(((BEI(i,:)*(xStep(i,end)-xStep(i,1))/100).^2 ))); %(xStep(i,end)-xStep(i,1))/100);
end

for i=1:1:length-1
    [r,lags] = xcorr(detrendCE(i,results.backwaterLengthIndex(i):shorelineIdx(i)),detrendCE(i+1,results.backwaterLengthIndex(i+1):shorelineIdx(i+1)));
    lagMax(i) = max(lags);
    rMaxCE(i) = max(xcorr(detrendCE(i,results.backwaterLengthIndex(i):shorelineIdx(i)),detrendCE(i+1,results.backwaterLengthIndex(i+1):shorelineIdx(i+1))));
    RMSce(i) = sqrt(sum(((detrendCE(i,results.backwaterLengthIndex(i):shorelineIdx(i))).^2 ))); %(xStep(i,end)-xStep(i,1))/100); *(xStep(i,end)-xStep(i,1))/100
end
rMaxCE(end) = rMaxCE(end-1);
rMaxCE = rMaxCE/max(rMaxCE);
lagMax(end) = lagMax(end-1);
lagMax = lagMax/max(lagMax);
RMSce = RMSce/max(RMSce);
fig(1) = figure(3);hold on;
set(gcf,'position',[100 100 300 250]);

p1 = plot(RMSce,'o','linewidth',1,'MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',0.1,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',2);
p2 = plot(rMaxCE,'o','linewidth',1,'MarkerFaceColor',[0.6 0.6 0.9],'LineWidth',0.1,'MarkerEdgeColor',[0.6 0.6 0.9],'MarkerSize',2);
p3 = plot(lagMax,'o','linewidth',1,'MarkerFaceColor',[0.9 0.6 0.6],'LineWidth',0.1,'MarkerEdgeColor',[0.9 0.6 0.6],'MarkerSize',2);
p4 = plot(movmean(RMSce,200),'k','linewidth',2);
p5 = plot(movmean(rMaxCE,200),'b','linewidth',2);
p6 = plot(movmean(lagMax,200),'r','linewidth',2);

plot([tAutobreak,tAutobreak],[0 1],':','color',[0.5 0.5 0.5],'lineWidth',1.5)
xlabel('time {\itt} (year)');
xlim([0 8000]);
xticks(0:2000:8000);
ylabel({'RMS,peak correlation','correlation lag'});
h1 = legend([p1 p2 p3 p4 p5 p6],{'RMS','PC','CL','averaged RMS','averaged PC','Averaged CL'},'location',...
    'southeast');

set(h1,'box','off');
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

scale = 0.8;
pos = get(gca, 'Position');
pos(4) = scale*pos(4);
set(gca, 'Position', pos)
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');
h = xlabel('{\itt*} ({\itt/T_a})');
Xtick = ((0:1:5)*T_auto);
Xticklabel = (num2str(0:1:5));
Xticklabel = split(Xticklabel);
xticks(Xtick);
xticklabels(Xticklabel);
xlim([0 8000]); %if not specify this, it will not plot properly
yticks([]);
yticklabels({});
ax2.XAxis.Color = 'k';
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

linkaxes([ax1,ax2])


