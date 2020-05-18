clear
clc

%Grid setup
noNodes = 200; % number of modeling nodes
dx = 5000; % x increment, unit is meter
dxp = dx; % previous dx for sed partition calculation
dxn = 1/noNodes; % distretized interval in moving boundary coordinates.
xn=0:dxn:1; %dimensionless grid

%simulation time-stepping
tf = 0.1; % ratio to adjust number of time steps and sea level rate
noTimeSteps = 8000/tf; % modeling time period
ti = 500/tf;  % sample interval for ploting figure
dt = 60*60*24*365*tf; % time increment of 1 year X tf

%sampling simulation results
noSkipSteps = 1/tf; %grab system information every noSkipSteps*dt time interval

%system thresholds (for autoretreat calculation)
rS_threshold = 0.1;
rS = ones(1,noNodes+1); %initial ratio of Sfrix/Sx

%channel system slope
Sfi = 2e-4; % initial slope of fluvial reach%system slopes
Sbb = 4.6e-4; % Slope of subaerial bedrock reach
Ssb = 1e-3; % slope of subaqueous basement
Sfore = 1.6e-3; % slope of foreset.
Sx = zeros(1,noNodes+1); % channel slope over x space

% fitted initial channel bed elevation and x node location
x = 0:dx:dx*noNodes;
x = x./1000;
xp = x; % x coordinates for previous time step
xShoreline = zeros(1,noNodes+1);% shoreline location
be = zeros(1,noNodes+1); %preallocate bed elevation
be(end) = -7; %bed elevation
be(1) = dx*noNodes*Sfi+be(end);%bed elevation
be = be(1):((be(end)-be(1))/noNodes):be(end);%bed elevation
bep = be; % bed elevation for previous time step
aggradationRate = be; % bed elevation increment/aggradation rate per year
bexc = ones(noTimeSteps,21);%preallocation for backwater length calculation

%inital basin configuration
seaLevel=0;
bebasei = -30; % initial bed elevation at foreset-basement break
sbai = -dx*noNodes; % initial x coordinate at bed rock/alluvial contact
stfi = 0; % initial x coordinate at shoreline break
ssbi = stfi+(be(end)-bebasei)/Sfore; % initial x coodinate of foreset bottom set boundary
sba = sbai; % x coordinate at bed rock/alluvial contact
stf = stfi; % x coordinate at shoreline break
ssb = ssbi; % x coordinate at foreset/bottom set contact.
bebase = be(end)-Sfore*(ssb-stf); % elevation of foreset/basement break, the value is updated each step, the value here is not important.
bebasep = bebase; %elevation of basement break for previous step, foreset increment calculation
autosb = nan(1,noTimeSteps); % index of shoreline break during shoreline break
autosbp = nan(1,noTimeSteps);
stfauto = []; % location of shoreline location during autoretreat

%basin rate variables
rsba = 0; % rate of migration of bedrock/alluvial contact 
rstf = 10; % rate of migration of shoreline break
rssb = 1; % rate of migration of foreset/bottom set contact

%storage variables
STF = zeros(1,noTimeSteps); % x coordinate at shoreline break throughout the modeling time frame
RSFT = zeros(1,noTimeSteps);% rate of migration of shoreline break for all the timed interval
STFauto = zeros(1,noTimeSteps); % collection of shoreline location during autoretreat

%preallocation of derivatives
dndt = zeros(1,noNodes+1); %elevation time der.
dqdxn = zeros(1,noNodes+1);  %sed flux time der.
dFdxn = zeros(11,noNodes+1); %grainsize classes time derv.
dLadxn = zeros(1,noNodes+1); %Active thickness of active layer time der
dDgdx = zeros(1,noNodes+1); %spatial change in median grainsize
Sx = zeros(1,noNodes+1); % slope over x space
Sfrix = zeros(1,noNodes+1); % friction slope over x space
dSdx = zeros(1,noNodes+1); % slope change over x space
dS = zeros(1,noNodes+1); % slope change over x space

x = xn.*(stf-sba)+sba;%initial node location

% hydraulics
un = 1.7;% uniform flow velocity m/s 
seaLevelRate = 10e-3*tf; % rate of sea level rise m/yr
If = 0.1;
Cf = 1/25^2; % friction coefficient
kc = 75e-3; % Roughness height
ar = 8.1; % Coefficient in Manning-Strickler resistance relation
g = 9.81; % gravitational acceleration
qw = (kc/ar^6*un^10/g^3/Sfi^3)^(1/4); % calculated qw
Qw = ones(1,noNodes+1).*qw;
bei = be; % set the initial bed elevation for comparison
sp = zeros(1,noNodes); % Flow surface profile
Hn = qw/un; % normal flow depth
qs = zeros(1,noNodes); % sediment flux
D50 = 120e-6; % medium grain size
R = 1.65; % Submerged specific gravity
au = 1; % explicit vs. implicit
LAMBDA = 1; % fraction of mud deposited in flood plain
Sinu = 1.7; % Sinuosity
rB = 60; % ratio between channel and flood plain
lambda = 0.4; % bed porosity

% flow depth calculation preallocation
H = zeros(1,noNodes+1); % flow depth.
H(end) = seaLevel-be(end); % flow depth this sets downstream boundary value as well
Hn = qw/un; % normal flow depth
H(1:noNodes) = Hn;
kh1 = zeros(1,noNodes+1);%4th oder runge-kutta method
kh2 = zeros(1,noNodes+1);
kh3 = zeros(1,noNodes+1);
kh4 = zeros(1,noNodes+1);
H_auto = zeros(1,noNodes+1);
H_auto_m = [];
t_auto_m = [];
det=0; % determine if it's time to set critical downstream flow depth

% grain size group and initial fractions for each group
GSGr = [-5 -4.5 -4 -3.5 -3 -2.5 -2 -1.5 -1 -0.5 0]; % grain size group
noGrainSizeGroups = length(GSGr);

%grain size calculation preallocation 
D = zeros(noGrainSizeGroups,noNodes+1); % grain size groups
DgPhi = zeros(1,noNodes+1); % geometric mean in phi scale
DgFsPhi = zeros(1,noNodes+1); % geometric mean in phi scale in substrate
Dv = zeros(1,noNodes+1); % Deviation of grain size
DvFs = zeros(1,noNodes+1); % Deviation of grain size in substrate
SDv = zeros(1,noNodes+1); % standard deviation
SDvFs = zeros(1,noNodes+1); % standard deviation
La = zeros(1,noNodes+1); % thickness of active layer
Lai = zeros(1,noNodes+1); % thickness of active layer for previous step
A = zeros(noGrainSizeGroups,noNodes+1); %
B = zeros(noGrainSizeGroups,noNodes+1); %
qit = zeros(1,noNodes+1); %
Chi = 1; % coefficient for calculating bedload layer sediment fraction
taog = ones(1,noNodes+1);
dSfdx = ones(1,noNodes+1);
dSxdx = ones(1,noNodes+1);

GSGr = [-5 -4.5 -4 -3.5 -3 -2.5 -2 -1.5 -1 -0.5 0]; % grain size group
GSGf = [0.6 1.7 4.4 9.2 15 38.2 15 9.2 4.4 1.7 0.6]/100; % grain size group fraction
Phi = zeros(11,noNodes+1); % Grain size group
F = zeros(11,noNodes+1); % fraction for each grain size group in active layer
f = zeros(11,noNodes+1); % fraction for each grain size group in bedload layer
for  ii=1:11
Phi(ii,:) = GSGr(ii);
F(ii,:) = GSGf(ii);
f(ii,:) = GSGf(ii);
end
Fs = F; % fraction each grain size group in substrate
Phia = sum(Phi(:,:).*F(:,:));
PhiaFs = sum(Phi(:,:).*Fs(:,:));
D(:,:) = 2.^Phi(:,:); % calculate mean grain size for each group in mm.
Dg = 2.^Phia/1000; % geometric mean grain size 
DgFs = 2.^PhiaFs/1000; % geometric mean grain size

% bws = 0; % backwater start point
% bwln = []; % backwater length using normal flow depth and slope
% bwlm = []; % backwater length using river mouth flow depth and slope
% bwlnx = []; % x coordinate of bw transition normal flow
% bwlmx = []; % x coordinate of bw transition river mouth

%pre-allocate for plotting
ti = 500/tf;  % sample interval for ploting figure

% pre-allocate model results output
%these fields are nxmxl; n = simulation number; m = grid nodes; 1 = time-step
noObservations = noTimeSteps/noSkipSteps;
tmp = nan(noObservations,noNodes+1);

flowDepth = tmp;
channelSlope = tmp;
aggradationRate = tmp;
grainSize = tmp;
sedimentFlux = tmp;
channelElevation = tmp;
nodeLocation = tmp;

backwaterLength = nan(1,noObservations);
backwaterLengthxIndex = nan(1,noObservations);
shorelineLocation = nan(1,noObservations);
shorelineLocation2 = nan(1,noObservations);
shorelineMigRate = nan(1,noObservations);
bwLengthNorFlow = nan(1,noObservations);%backwater length based on slope and flow depth at normal flow
bwLengthMouth = nan(1,noObservations);%backwater length based on slope and flow depth at river mouth
bwLengthSL = nan(1,noObservations);%backwater length channel bed intercepts sea level
sedPartition =  nan(1,noObservations); 
basementPosition = nan(2,noObservations); 
maxAggradationX = nan(1,noObservations); 
shorelineIdx = nan(1,noObservations);
ic = 1;
det = 1;

% backwater calculation
for tspan=1:noTimeSteps
    idxu = (1:noNodes+1); % node index upstream of auto shoreline
    idxd = (1:noNodes+1); % node index downstream of auto shoreline
    idxu = idxu(rS>rS_threshold);
    idxd = idxd(rS<=rS_threshold);
    [colum, idxul]= size(rS(rS>rS_threshold));
    idx = (1:noNodes);
    
    kh1(idx) = ((-be(idx+1)+be(idx))/dxn/(stf-sba)...
                -(kc^(1/3)*qw^2/ar^2/g./H(idx+1).^(10/3)))./(1-qw^2/g./H(idx+1).^3)*dxn*(stf-sba);      
    kh2(idx) = ((-be(idx+1)+be(idx))/dxn/(stf-sba)...
                -(kc^(1/3)*qw^2/ar^2/g./(H(idx+1)-kh1(idx)/2).^(10/3)))./(1-qw^2/g./(H(idx+1)-kh1(idx)/2).^3)*dxn*(stf-sba);        
    kh3(idx) = ((-be(idx+1)+be(idx))/dxn/(stf-sba)...
                -(kc^(1/3)*qw^2/ar^2/g./(H(idx+1)-kh2(idx)/2).^(10/3)))./(1-qw^2/g./(H(idx+1)-kh2(idx)/2).^3)*dxn*(stf-sba);  
    kh4(idx) = ((-be(idx+1)+be(idx))/dxn/(stf-sba)...
                -(kc^(1/3)*qw^2/ar^2/g./(H(idx+1)-kh3(idx)).^(10/3)))./(1-qw^2/g./(H(idx+1)-kh3(idx)).^3)*dxn*(stf-sba);  
    H(idx) = H(idx+1)-kh1(idx)/6-kh2(idx)/3-kh3(idx)/3-kh4(idx)/6;    
    H(H<=Hn) = Hn; % upstream flow depth might get below Hn due to model instability        
    qs = (R*g*Dg.^3).^0.5./(R*g*Dg).^3*0.0355/ar^4.*(kc./H).^(2/3).*(qw./H).^6; % calculate qs
    
    autosb(tspan) = idxul;%storing location index of shoreline
    if idxu(end)~=noNodes+1
        H_auto(tspan) = interp1(rS(autosb(tspan):end),H(autosb(tspan):end),rS_threshold,'spline');
        if det == 1
            tAutobreak = tspan;%determine the timing of autobreak
            det = 0;
        end
    else
        H_auto(tspan) = H(idxu(end));%storing flow depth at shoreline break
    end
    Qw = qw;
   
    % grain size calculation
    taog = 1/(ar^2)*(kc./H).^(1/3)*qw^2/R/g./Dg./H.^2;
    % calculate coefficient A and B from Naito et al. (2019)
    A(:,:) = 0.455*(D(:,:)./Dg).^(-0.839);
    B(:,:) = 0.353*(D(:,:)./Dg).^(-1.157);
    qit = sum(F(:,:).*A(:,:).*(taog.*Dg./D(:,:)).^B(:,:));

    % calculate fraction of grain size in bedload
    f(:,:) = F(:,:).*A(:,:).*(taog.*Dg./D(:,:)).^B(:,:)./qit;    
    La = 2*8*H.*(Dg./H).^0.3; %Dg(M) * 2^(1.28*SDv(M)); % thickness of active layer
    % calculate fraction in substrate
    Fs(:,:) = Chi * F(:,:) + (1-Chi) * f(:,:);

    % end of grain size calculation        
    stfauto = x(autosb(tspan));
    STFauto(tspan) = stfauto;

    if tspan==1
        dqdxn(end) = (qs(end)-qs(end-1))/dxn;
        dqdxn(1) = (qs(2)-qs(1))/dxn;
        dndt(end) = 1e-5*(-Sfi)-dqdxn(end)*If*(1+LAMBDA)*Sinu/rB/(1-lambda)/(stf-sba);
        dndt(1) = 1e-5*(-Sfi)-dqdxn(1)*If*(1+LAMBDA)*Sinu/rB/(1-lambda)/(stf-sba);
    end

    rstf = 1/Sfore*(If*Sinu*(1+LAMBDA)*qs(end)/rB/(1-lambda)/(ssb-stf)-dndt(end)); % rate of migration of shoreline break
    rssb = 1/(Sfore-Ssb)*(Sfore*rstf+dndt(end)); % rate of rate of migration of foreset/bottom set contact
    rsba = -1/Sbb*dndt(1); % rate of migration of bedrock/alluvial contact

    if qs(end)==0
        rstf = 0;
        rssb = 0;
    end
    
    RSFT(tspan) = rstf;
    Sxp = Sx(end); % Sx at downstream end from previous iteration, for sediment partition calculation 
    bet = be(end); % temporary be(end) for mass balance calculation
    Sxendi = Sx(end); % downstream slope at previous time step

    %calculate spatial derivatives of sediment flux, grain size fraction,
    %flow depth and active layer thickness
    %downstream point
    dqdxn(noNodes+1) = (qs(noNodes+1)-qs(noNodes))/dxn;
    dndt(noNodes+1) = (xn(noNodes+1)*rstf+(1-xn(noNodes+1))*rsba)/(stf-sba)*(be(noNodes+1)-be(noNodes))/dxn...
    -dqdxn(noNodes+1)*If*(1+LAMBDA)*Sinu/rB/(1-lambda)/(stf-sba);
    dFdxn(:,noNodes+1) = (F(:,noNodes+1)-F(:,noNodes))/dxn;
    dLadxn(noNodes+1) = (La(noNodes+1)-La(noNodes))/dxn;
    dHdx(noNodes+1) = (H(noNodes+1)-H(noNodes))/(x(noNodes+1)-x(noNodes));

    % upstream point
    dqdxn(1) = (qs(2)-qs(1))/dxn;
    dndt(1) = (xn(1)*rstf+(1-xn(1))*rsba)/(stf-sba)*(be(2)-be(1))/dxn...
    -dqdxn(1)*If*(1+LAMBDA)*Sinu/rB/(1-lambda)/(stf-sba);
    dFdxn(:,1) = (F(:,2)-F(:,1))/dxn;%dFdxn(ii,M) = (F(ii,M)-GSGf(ii))/dxn;
    dLadxn(1) = (La(2)-La(1))/dxn;
    dHdx(1) = (H(2)-H(1))/(x(2)-x(1));

    % cases 2 to N points
    dqdxn(2:noNodes) = au*(qs(2:noNodes)-qs(1:noNodes-1))/dxn+(1-au)*(qs(3:noNodes+1)-qs(1:noNodes-1))/2/dxn;
    dndt(2:noNodes) = (xn(2:noNodes)*rstf+(1-xn(2:noNodes))*rsba)/(stf-sba).*(be(2:noNodes)-be(1:noNodes-1))/dxn...
    -dqdxn(2:noNodes)*If*(1+LAMBDA)*Sinu/rB/(1-lambda)/(stf-sba);
    dFdxn(:,2:noNodes) = au*(F(:,2:noNodes)-F(:,1:noNodes-1))/dxn+(1-au)*(F(:,3:noNodes+1)-F(:,1:noNodes-1))/2/dxn;
    dLadxn(2:noNodes) = au*(La(2:noNodes)-La(1:noNodes-1))/dxn+(1-au)*(La(3:noNodes+1)-La(1:noNodes-1))/2/dxn;
    dHdx(2:noNodes) = (H(2:noNodes)-H(1:noNodes-1))./(x(2:noNodes)-x(1:noNodes-1));
    
    %update bed elevation and flow surface profile
    be = be+dndt*dt;
    sp = be+H;

    Sx(noNodes+1) = -(be(noNodes+1)-be(noNodes))/(x(noNodes+1)-x(noNodes));
    Sx(1) = -(be(2)-be(1))/(x(2)-x(1));
    Sx(2:noNodes) = -(be(2:noNodes)-be(1:noNodes-1))./(x(2:noNodes)-x(1:noNodes-1));
    Sfrix = kc^(1/3)/ar^2.*Qw.^2./g./H.^(10/3);

    dSfdx(noNodes+1) = (Sfrix(noNodes+1)-Sfrix(noNodes))/(x(noNodes+1)-x(noNodes));
    dSfdx(1) = (Sfrix(2)-Sfrix(1))/(x(2)-x(1));
    dSfdx(2:noNodes) = (Sfrix(2:noNodes)-Sfrix(1:noNodes-1))./(x(2:noNodes)-x(1:noNodes-1));
    dSxdx(noNodes+1) = (Sx(noNodes+1)-Sx(noNodes))/(x(noNodes+1)-x(noNodes));
    dSxdx(1) = (Sx(2)-Sx(1))/(x(2)-x(1));
    dSxdx(2:noNodes) = (Sx(2:noNodes)-Sx(1:noNodes-1))./(x(2:noNodes)-x(1:noNodes-1));
    rS = Sfrix./Sx;
    
    % calculate updated fractions of grain size in active layer     
    idx = (1:noNodes+1);
    idx = idx(qs~=0); % no calculation for drowned part (qs=0) of the channel
    F(:,idx) = F(:,idx)+(xn(idx)*rstf+(1-xn(idx))*rsba)/(stf-sba).*dFdxn(idx)*dt...
    -1./La(idx).*(F(:,idx)-Fs(:,idx)).*((La(idx)-Lai(idx))/dt...
    -(xn(idx)*rstf+(1-xn(idx))*rsba)/(stf-sba).*dLadxn(idx))*dt...
    -If*(1+LAMBDA)*Sinu./(La(idx)*(1-lambda)*(stf-sba)*rB).*(dqdxn(idx).*f(:,idx)-Fs(:,idx).*dqdxn(idx))*dt;

    F(F<0) = eps;
    Lai=La;% update old active layer Lai to new La for next cycle
    Dv=zeros(noNodes+1); % reset the deviation to be zero for calculation of next node.

    Dg = sum(D(:,:).*F(:,:))/1000; % geometric mean grain size 
    DgFs = sum(D(:,:).*Fs(:,:))/1000; % geometric mean grain size 

    seaLevel = seaLevel+seaLevelRate;
    H(end)=seaLevel-be(end);% update river mouth flow depth

    x = xn.*(stf-sba)+sba; % convert normalized x distance back to dimensional
    
    ssbp = ssb; % temporary ssb for mass balance calculation
    stft = stf; % temporary stf for mass balance calculation   
    ssbt = ssb; % temporary ssb for mass balance calculation
    stft = stf; % temporary stf for mass balance calculation
    bebaset = bebase; % temporary bebase for mass balance calculation

    sba = sba+rsba*dt; % update sba
    ssb = ssb+rssb*dt; % update ssb
    stf = stf+rstf*dt; % update stf   
    STF(tspan) = stf;
    ssby = be(end)-Sfore*(ssb-stf);%elevation of basement foreset break

    % calculate bed elevation increment
    if rstf>0
        BEInt = interp1(x,be,xp,'spline');  
        BEI = BEInt-bep;
    else
        BEInt = interp1(xp,bep,x,'spline');
        BEI = be - BEInt;
    end
    [depoFront,idx] = max(BEI);
    depoFrontX = x(idx);
    bep = be; % update BEP for next round of iteration

    % Calculate location of BW transition
    idx = (100:noNodes+1);
    idx = idx(dHdx(100:noNodes+1)<5e-6);
    bwl(tspan) = x(autosb(tspan))-x(idx(end));
    bwlxIdx = idx(end);
    bexc(tspan,:) = interp1(x,be,x(autosb(tspan))-200000:10000:x(autosb(tspan)),'spline');

    
    % Mass Balance calculation
    if rstf>0
        fa1 = sum(BEI) * dxp - 0.5*BEI(end)*dxp;
        if be(end)<bep(end)
            fa2 = 0.5 * Sx(end) * (x(end)-xp(end))^2 - 0.5 * (bep(end)-be(end))*(x(end)-xp(end));
            da1 = 0.5 * (x(end)-xp(end)-(bep(end)-be(end))/Sfore) * (bep(end)-be(end));
            da2 = 0.5 * (ssb-ssbp - (bebasep-bebase)/Sfore) * (bep(end)-bebasep);
        else
            fa2 = 0.5 * Sx(end) * (x(end)-xp(end))^2 + 0.5 * (be(end)-bep(end)) * (x(end)-xp(end));
            da1 = 0.5 * (x(end)-xp(end)+(bep(end)-be(end))/Sfore) * (be(end)-bep(end));
            da2 = 0.5 * (ssb-ssbp - (bebasep-bebase)/Sfore) * (bep(end)-bebasep);
        end
        else  %rstf<0 shoreline retreat
            fa1 = sum(BEI) * dx - 0.5*BEI(end)*dx;
            fa2 = 0.5 * (xp(end)-x(end)) * (be(end)-bep(end)) - 0.5 * Sxp * (xp(end)-x(end))^2;
            da1 = 0.5 * (ssb-ssbp - (bebasep-bebase)/Sfore) * (be(end)-bep(end));
            da2 = 0.5 * (ssb-ssbp - (bebasep-bebase)/Sfore) * (bep(end)-bebasep);  
    end
    da3 = 0.5 * (bebasep-bebase) * ((ssb-ssbp) - (bebasep-bebase)/Sfore);
    fa = fa1 + fa2;
    da = da1 + da2 + da3;
    fp = fa / (fa+da);
    sedPar = fp;
    
    xp = x; % update xp so for next iteration of calculation, this is placed after ploting because plotting needs previous x location.
    bep = be; % update BEP so for next iteration of calculation
    bebasep = bebase; % temporary bebase for mass balance calculation

    % plotting
    if tspan==1||tspan/500*tf-round(tspan/500*tf)==0
    figure(1);hold on,box on
    hold on
    aa=plot(x./1e3,be,'k');
    p1 = plot(stf/1e3,be(end),'bo');

    if rssb>0
    cc=plot([x(end)/1e3,ssb/1e3],[be(end),ssby],'k'); % plot delta foreset
    end
    xlabel('distance (km)', 'FontSize', 10);
    ylabel('elevation (m)','FontSize', 10);

    if rS(end)<rS_threshold
    p2 = plot(stfauto/1e3,be(autosb(tspan)),'ro');
    legend([p1,p2],{'shoreline before autobreak','shoreline after autobreak'});
    end
    xlim([-400 100])
    
    figure(2) % grain size in bedload material load 
    hold on
    plot(xn,Dg)
    xlabel('distance (m)', 'FontSize', 10);
    ylabel('grain size in bedload (m)','FontSize', 10);

    figure(3)
    subplot(1,2,1);hold on
    plot(xn,H,'b')
    xlabel('distance (m)', 'FontSize', 10);ylabel('flow depth (m)','FontSize', 10);

    subplot(1,2,2);hold on
    plot(x,rS,'o')
    xlabel('distance (m)', 'FontSize', 10);ylabel('rS','FontSize', 10);
    end
    xp = x; % update xp so for next iteration of calculation, this is placed after ploting because plotting needs previous x location.
    
    %store results
    if mod(tspan,noSkipSteps) == 0
        flowDepth(ic,:) = H;
        channelSlope(ic,:) = Sx;
        aggradationRate(ic,:) = BEI;
        grainSize(ic,:) = Dg;
        sedimentFlux(ic,:) = qs;
        channelElevation(ic,:) = be;
        nodeLocation(ic,:) = x;
        
        backwaterLength(ic) = bwl(tspan);
        backwaterLengthxIndex(ic) = bwlxIdx;
        shorelineLocation(ic) = stf;
        shorelineLocation2(ic) = stfauto;
        shorelineIdx(ic) = autosb(tspan);
        shorelineMigRate(ic) = rstf;
        sedPartition(ic) = sedPar;
        basementPosition(1,ic) = ssb; %x coordinate of basement foreset break
        basementPosition(2,ic) = ssby; %y coordinate of basement foreset break
        maxAggradationX(ic) = depoFrontX; % location of maximum deposition
        ic = ic + 1; 
    end

end
tAutobreak = tAutobreak*tf;
shorelineMigRate2 = (shorelineLocation2(2:end)-shorelineLocation2(1:end-1));
shorelineMigRate2 = [shorelineMigRate2(1),shorelineMigRate2];

%return values in the form of a structure
results.flowDepth = flowDepth;
results.channelSlope = channelSlope;
results.aggradationRate = aggradationRate;
results.maxAggradationX = maxAggradationX;
results.grainSize = grainSize;
results.sedimentFlux = sedimentFlux;
results.channelElevation = channelElevation;
results.nodeLocation = nodeLocation;

results.backwaterLength =  backwaterLength;
results.backwaterLengthIndex =  backwaterLengthxIndex;
results.shorelineLocation =  shorelineLocation;
results.shorelineLocation2 =  shorelineLocation2;
results.shorelineMigRate = shorelineMigRate;
results.sedPartition = sedPartition;
results.basementPosition = basementPosition;
results.tAutobreak = tAutobreak;
results.shorelineIdx = shorelineIdx;

shorelineInterp = shorelineLocation2(tAutobreak:end);
shorelineInterpT = (tAutobreak:5000);
shorelineInterp = -0.002519*shorelineInterpT.^2 - 25.62*shorelineInterpT + 6.252e4;
shorelineInterp = [shorelineLocation2(1:tAutobreak-1),shorelineInterp ];
shorelineMigRate2 = (shorelineInterp(2:end)-shorelineInterp(1:end-1));
shorelineMigRate2 = [shorelineMigRate2(1),shorelineMigRate2];
shorelineMigRate2 (tAutobreak) = shorelineMigRate2 (tAutobreak+1);
dMigRate = shorelineMigRate2(tAutobreak+1:end);
dMigRate = dMigRate(2:end)-dMigRate(1:end-1);

idx = (1:5000);
idx = idx(shorelineMigRate2>0);
tRetreat = idx(end);%time retreat starts

%calculate autostratigraphic time and length scales
L_auto = qs(1)/10e-3*60*60*24*365/60/10*2*1.7/1e3/0.6;%autostratiraphic length scale
T_auto = Sfi * qs(1)/(10e-3)^2/60/10*2*1.7/0.6*(60*60*24*365);%autostratigraphic time scale

% plotting results
figure(4)
hold on
yyaxis left
plot(x,H);xlabel('distance (m)', 'FontSize', 10);ylabel('flow depth (m)','FontSize', 10);
yyaxis right
plot(x,Dg*1e6);ylabel('grain size (\mum)','FontSize', 10);

figure(1)
hold on
plot(x./1e3,be+H,'-b')
plot([stf/1e3,ssb/1e3],[seaLevel,seaLevel],'-b')

close all
