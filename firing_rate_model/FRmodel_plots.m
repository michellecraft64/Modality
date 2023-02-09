% script for firing rate model in paper
% has 1 excit, 1 inhib; implem in WCnet_EIcells.m (or WCnet_EIcellsF.m)

%% 
clear

flag_Regime=1; %1=Syn Depr;  2=Facilitation

%Select params
tdE_arr=10; %in ms
trE_arr=2; 
tw=100; %in ms
w_Parr=[0.2 0.25 0.15]; %0.2 original, 0.25 MUS, 0.15 BIC
crr=1; % crr param; turn off plasticity with crr=0
plot_Flag=1; %Flag to plot comparisons of parameter differences
NBdAcc_Flag=0; %Flag to sim Negative Binomial trial to trial differences

ParamMat=combvec([tdE_arr;trE_arr],w_Parr);
numParams=size(ParamMat,2);

%Time variables
TimeVars=struct('dt',0.01,'lenSpon',0,'lenEvok',900,'lenTrn',500); %opt Twin = 900ms evoked (2s spont)

% % weights below & lamO/lamR 
%with original lamO/lamR
alphMaxRet_MC=200; %Constant multiplier for input lam_Ret to MC
alphMaxOr_MC=100; %Constant multiplier for input lam_Or to MC
alphMaxRet_I=120; %Constant multiplier for input lam_Ret to I
alphMaxOr_I=150; %Constant multiplier for input lam_Or to I

%Define time vectors
dt=TimeVars.dt;
tFullEnd=TimeVars.lenSpon+TimeVars.lenEvok+TimeVars.lenTrn;
tEnd=TimeVars.lenSpon+TimeVars.lenEvok; %without transient time
t=(0:dt:tEnd)';
Lt=length(t); %same length as Tmev
trn=(0:dt:tFullEnd-tEnd-dt)';
lenTrn=length(trn); %transient time
tmeEvok=(0:dt:TimeVars.lenEvok)';
tmeSpon=(-TimeVars.lenSpon:dt:0)';
Tmev=[tmeSpon(1:end-1); tmeEvok]/1000;
lenFullTme=tFullEnd/dt+1;
%Ortho input
lamO=struct('tsft',25.1,'lmOsp',.002,'lmOevk',15/367,'tauO',60,'tauO2',90,'ssV',.007); % derived from Craft et al 21
tmShift=50; %in ms
TimeVars.lenEvokS=[tmeEvok; tmeEvok(end)+(dt:dt:tmShift)']; %take nuevok out further by tmShift
nmShift=length((dt:dt:tmShift)); %# elem must shift by
nuevokO=lamO.lmOevk*(-(TimeVars.lenEvokS+lamO.tsft).*exp(-(TimeVars.lenEvokS+lamO.tsft)/lamO.tauO) + (TimeVars.lenEvokS+lamO.tsft).*exp(-(TimeVars.lenEvokS+lamO.tsft)/lamO.tauO2))+lamO.ssV;
nuevokO=nuevokO(nmShift+1:end); %shift over by tmShift
lam_O=[lamO.lmOsp*ones(length(tmeSpon)-1,1);nuevokO];
%Retro input
lamR=struct('tsft',25.1,'lmRspon',.002,'lmRevk',0.3/367,'tau_R',450); % derived from Craft et al 21
nuevok_R=(tmeEvok+lamR.tsft).*exp(-(tmeEvok+lamR.tsft)/lamR.tau_R);
lam_R=[lamR.lmRspon*ones(length(tmeSpon)-1,1);lamR.lmRevk*nuevok_R];

InputCurrent_Or=struct('Mc',lam_O*alphMaxOr_MC,'Pg',lam_O*alphMaxOr_I);
InputCurrent_Ret=struct('Mc',lam_R*alphMaxRet_MC,'Pg',lam_R*alphMaxRet_I);

xMin=struct('xM',0,'xP',0); %minimal current for firing

%Initialize output
Am_O=zeros(numParams,Lt);
Ap_O=zeros(numParams,Lt);
syEo=zeros(numParams,Lt);
syIo=zeros(numParams,Lt);
w_MP_O=zeros(numParams,Lt);
Am_R=zeros(numParams,Lt);
Ap_R=zeros(numParams,Lt);
syEr=zeros(numParams,Lt);
syIr=zeros(numParams,Lt);
w_MP_R=zeros(numParams,Lt);
% save transients too
Amt_O=zeros(numParams,lenTrn);
Apt_O=zeros(numParams,lenTrn);
w_MPt_O=zeros(numParams,lenTrn);
Amt_R=zeros(numParams,lenTrn);
Apt_R=zeros(numParams,lenTrn);
w_MPt_R=zeros(numParams,lenTrn);
netSpikRt0=zeros(numParams,2);


%Iterate through params
for i=1:numParams
    %Define specific time constants and coupling strengths
    TimeConstants=struct('tau_M',10,'tau_P',5.5,'tdE',ParamMat(1,i),'trE',ParamMat(2,i),'tdI',10,'trI',2,'tw',tw);
    % td[E/I]=decay time of syn, tr[E/I]=rise time of syn
    CpldWx=struct('w_I',1,'w_P',ParamMat(3,i),'wmp0',1); 
    
    switch flag_Regime
        case 1
            %call for Ortho
            [Am_O(i,:),Ap_O(i,:),syEo(i,:),syIo(i,:),w_MP_O(i,:),Amt_O(i,:),Apt_O(i,:),w_MPt_O(i,:)]=WCnet_EIcells(CpldWx,TimeVars,TimeConstants,xMin,InputCurrent_Or,crr);
            %call for Ret
            [Am_R(i,:),Ap_R(i,:),syEr(i,:),syIr(i,:),w_MP_R(i,:),Amt_R(i,:),Apt_R(i,:),w_MPt_R(i,:)]=WCnet_EIcells(CpldWx,TimeVars,TimeConstants,xMin,InputCurrent_Ret,crr);
        case 2
            %call for Ortho
            [Am_O(i,:),Ap_O(i,:),syEo(i,:),syIo(i,:),w_MP_O(i,:),Amt_O(i,:),Apt_O(i,:),w_MPt_O(i,:)]=WCnet_EIcellsF(CpldWx,TimeVars,TimeConstants,xMin,InputCurrent_Or,crr);
            %call for Ret
            [Am_R(i,:),Ap_R(i,:),syEr(i,:),syIr(i,:),w_MP_R(i,:),Amt_R(i,:),Apt_R(i,:),w_MPt_R(i,:)]=WCnet_EIcellsF(CpldWx,TimeVars,TimeConstants,xMin,InputCurrent_Ret,crr);
    end

    %net spike rate
    netSpikRt0(i,:)=[sum(Am_O(i,:))*dt sum(Am_R(i,:))*dt]./(1000*Tmev(end)); %b/c Tmev in s & nu0 is per ms, convt to Hz!
end



%% Neg bin decode acc
aspPv=[0.12 0.12; .12 .04; .12 .12];
if NBdAcc_Flag==1
    tic
    Ntr=50000; %aspirational p for neg binomial, p=1 is Poisson, p inverse prop to variance
    AccMat=zeros(numParams,4);
    for k=1:numParams
        
        nuO=squeeze(Am_O(k,:));
        % Negative Binomial
        meanOc=sum(nuO)*dt/(1000*Tmev(end)); %b/c Tmev in s & nu0 is per ms, convt to Hz!
        sc(k,1)=aspPv(k,1)/(1-aspPv(k,1)); %r=sc*mu, sc=p/(1-p) in neg binomials
        rO=meanOc*sc(k,1);
        
        nuR=squeeze(Am_R(k,:));
        meanRc=sum(nuR)*dt/(1000*Tmev(end)); %b/c Tmev in s & nu0 is per msm convt to Hz!
        sc(k,2)=aspPv(k,2)/(1-aspPv(k,2));
        rR=meanRc*sc(k,2);
        sumSpks_O_Close=nbinrnd(rO,rO/(meanOc+rO),Ntr,1);
        sumSpks_R_Close=nbinrnd(rR,rR/(meanRc+rR),Ntr,1);

        figure('Renderer', 'Painters');
        hold on
        histogram(sumSpks_O_Close,'Normalization','probability','FaceColor','b','EdgeColor','b','LineStyle','none')
        histogram(sumSpks_R_Close,'Normalization','probability','FaceColor','r','EdgeColor','r','LineStyle','none')
        box off
        %Set decode acc params
        allRespn=[sumSpks_O_Close;sumSpks_R_Close];
        mnOr=mean(sumSpks_O_Close,1);
        mnRet=mean(sumSpks_R_Close,1);
        %Output (trial-avg net responses (netRespn); all LDA accuracies (OrAcc/RetAcc); averaged LDA accuracies (NeuronClass))
        % ---- speed up code by only checking out sub-set of threshold ---
        if(mnRet>mnOr) 
            lwThr=prctile(sumSpks_O_Close,35); %35th percentile of smaller
            hiThr=prctile(sumSpks_R_Close,65); %65th percentile of larger
        else
            lwThr=prctile(sumSpks_R_Close,35); %35th percentile of smaller
            hiThr=prctile(sumSpks_O_Close,65); %65th percentile of larger
        end
        indT=(allRespn>=lwThr)&(allRespn<=hiThr);
        thrsh=allRespn(indT); %consider subset for threshold values
        Nthrs=length(thrsh);
        RetAcc=zeros(Nthrs,1);
        OrAcc=zeros(Nthrs,1);
        ttlAcc=zeros(Nthrs,1);
        for i=1:Nthrs %consider subset for threshold values
            Respn=ones(Ntr*2,1);
            if mnOr>mnRet
                Respn(allRespn>thrsh(i))=0;
            elseif mnRet>mnOr
                Respn(allRespn<thrsh(i))=0;
            else
                Respn(allRespn<thrsh(i))=0;
            end
            OrAcc(i)=1-sum(Respn(1:Ntr))/Ntr;
            RetAcc(i)=sum(Respn(Ntr+1:end))/Ntr;
            ttlAcc(i)=(OrAcc(i)+RetAcc(i))/2;
        end % ---- end of 'new code' ----
        [~,OptInd]=max(ttlAcc);
        AccMat(k,:)=[thrsh(OptInd) OrAcc(OptInd) RetAcc(OptInd) ttlAcc(OptInd)];
        % plot optim threshold
        my=gca;
        my=my.YLim(end);
        Yvc=(0: (my)/(100-1) : my)';
        plot(thrsh(OptInd)*ones(100,1),Yvc,'k','LineWidth',2)
        
    end %over all params
    AccMat(:,4) %show all decoding accuracies (total)
    toc
end

%% plots %%
if(plot_Flag)
%Plot different "drugs"
%MC - Ortho
ccD=[128.01 128.01 128.01; 140 0 255; 0 200 0]./255;
figure('Renderer', 'Painters')
hold on
plot(Tmev,Am_O(1,:),'LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Am_O(3,:),'LineWidth',2,'Color',ccD(3,:))
plot(Tmev,Am_O(2,:),'LineWidth',2,'Color',ccD(2,:))
legend('No Drug','Bicuculine','Muscimol')
yEnd=ylim;
ylabel('MC Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'FontSize',18)
set(gca,'XLim',[Tmev(1) Tmev(end)])
load dBase_noDepr.mat
plot(Tmev,Am_base_O(1,:),'--','LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Am_base_O(2,:),'--','LineWidth',2,'Color',ccD(2,:))
plot(Tmev,Am_base_O(3,:),'--','LineWidth',2,'Color',ccD(3,:))

%MC - Retro
figure('Renderer', 'Painters')
hold on
plot(Tmev,Am_R(1,:),'LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Am_R(3,:),'LineWidth',2,'Color',ccD(3,:))
plot(Tmev,Am_R(2,:),'LineWidth',2,'Color',ccD(2,:))
legend('No Drug','Bicuculine','Muscimol')
ylim(yEnd)
ylabel('MC Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'FontSize',18)
set(gca,'XLim',[Tmev(1) Tmev(end)])
plot(Tmev,Am_base_R(1,:),'--','LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Am_base_R(2,:),'--','LineWidth',2,'Color',ccD(2,:))
plot(Tmev,Am_base_R(3,:),'--','LineWidth',2,'Color',ccD(3,:))


%I - Ortho
figure('Renderer', 'Painters')
hold on
plot(Tmev,Ap_O(1,:),'LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Ap_O(3,:),'LineWidth',2,'Color',ccD(3,:))
plot(Tmev,Ap_O(2,:),'LineWidth',2,'Color',ccD(2,:))
legend('No Drug','Bicuculine','Muscimol')
yEnd=ylim;
ylabel('PGC Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'FontSize',18)
set(gca,'XLim',[Tmev(1) Tmev(end)])
plot(Tmev,Ap_base_O(1,:),'--','LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Ap_base_O(2,:),'--','LineWidth',2,'Color',ccD(2,:))
plot(Tmev,Ap_base_O(3,:),'--','LineWidth',2,'Color',ccD(3,:))

%I - Retro
figure('Renderer', 'Painters')
hold on
plot(Tmev,Ap_R(1,:),'LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Ap_R(3,:),'LineWidth',2,'Color',ccD(3,:))
plot(Tmev,Ap_R(2,:),'LineWidth',2,'Color',ccD(2,:))
legend('No Drug','Bicuculine','Muscimol')
ylim(yEnd)
ylabel('I Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'FontSize',18)
set(gca,'XLim',[Tmev(1) Tmev(end)])
plot(Tmev,Ap_base_R(1,:),'--','LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Ap_base_R(2,:),'--','LineWidth',2,'Color',ccD(2,:))
plot(Tmev,Ap_base_R(3,:),'--','LineWidth',2,'Color',ccD(3,:))

%Plot weights and thresholds over time
%Ortho
figure('Renderer', 'Painters')
hold on
for j=1:3
    plot(Tmev,CpldWx.wmp0+w_MP_O(j,:),'color',ccD(j,:),'LineWidth',2)
end
plot(Tmev,CpldWx.wmp0*ones(size(Tmev)),'k--','LineWidth',2)
xlabel('Time (s)')
ylabel('Syn Weights','FontSize',12)
title('Ortho','FontSize',18)
set(gca,'XLim',[Tmev(1) Tmev(end)])
%Retro
figure('Renderer', 'Painters')
hold on
for j=1:3
    plot(Tmev,CpldWx.wmp0+w_MP_R(j,:),'color',ccD(j,:),'LineWidth',2)
end
plot(Tmev,CpldWx.wmp0*ones(size(Tmev)),'k--','LineWidth',2)
xlabel('Time (s)')
ylabel('Syn Weights','FontSize',12)
set(gca,'XLim',[Tmev(1) Tmev(end)])
title('Retro','FontSize',18)

%table of relative change
or_nd_Dep=sum(Am_O(1,:))*dt/(1000*Tmev(end));
or_bic_Dep=sum(Am_O(3,:))*dt/(1000*Tmev(end));
or_mus_Dep=sum(Am_O(2,:))*dt/(1000*Tmev(end));
    rt_nd_Dep=sum(Am_R(1,:))*dt/(1000*Tmev(end));
    rt_bic_Dep=sum(Am_R(3,:))*dt/(1000*Tmev(end));
    rt_mus_Dep=sum(Am_R(2,:))*dt/(1000*Tmev(end));
or_nd_base=sum(Am_base_O(1,:))*dt/(1000*Tmev(end));
or_bic_base=sum(Am_base_O(3,:))*dt/(1000*Tmev(end));
or_mus_base=sum(Am_base_O(2,:))*dt/(1000*Tmev(end));
    rt_nd_base=sum(Am_base_R(1,:))*dt/(1000*Tmev(end));
    rt_bic_base=sum(Am_base_R(3,:))*dt/(1000*Tmev(end));
    rt_mus_base=sum(Am_base_R(2,:))*dt/(1000*Tmev(end));

drgComp={'ND/Bic';'ND/Mus'}; 
orth_NoDeprs=[or_bic_base-or_nd_base ; or_nd_base-or_mus_base  ]./or_nd_base;
orth_Deprs=[or_bic_Dep-or_nd_Dep ; or_nd_Dep-or_mus_Dep  ]./or_nd_Dep;
ret_NoDeprs=[rt_bic_base-rt_nd_base ; rt_nd_base-rt_mus_base  ]./rt_nd_base;
ret_Deprs=[rt_bic_Dep-rt_nd_Dep ; rt_nd_Dep-rt_mus_Dep  ]./rt_nd_Dep;
T=table(drgComp,orth_NoDeprs,orth_Deprs,ret_NoDeprs,ret_Deprs)

end

if (0)
    %% Transient plots .  Change if(0) to if(1) if you want to see, or just run this piece of code
    %Plot different "drugs"
    %MC - Ortho
    figure('Renderer', 'Painters')
    hold on
    plot(trn/1000,Amt_O(1,:),'LineWidth',2,'Color',[160/255 160/255 160/255])
    plot(trn/1000,Amt_O(3,:),'LineWidth',2,'Color',[10/255 206/255 10/255])
    plot(trn/1000,Amt_O(2,:),'LineWidth',2,'Color',[144/255 0 253/255])
    legend('No Drug','Bicuculine','Muscimol')
    yEnd=ylim;
    ylabel('MC Firing Rate (Hz)')
    xlabel('Time (s)')
    title('Ortho')
    set(gca,'FontSize',18)
    set(gca,'XLim',[trn(1) trn(end)]./1000)
    %MC - Retro
    figure('Renderer', 'Painters')
    hold on
    plot(trn/1000,Amt_R(1,:),'LineWidth',2,'Color',[160/255 160/255 160/255])
    plot(trn/1000,Amt_R(3,:),'LineWidth',2,'Color',[10/255 206/255 10/255])
    plot(trn/1000,Amt_R(2,:),'LineWidth',2,'Color',[144/255 0 253/255])
    legend('No Drug','Bicuculine','Muscimol')
    ylim(yEnd)
    ylabel('MC Firing Rate (Hz)')
    xlabel('Time (s)')
    title('Retro')
    set(gca,'FontSize',18)
    set(gca,'XLim',[trn(1) trn(end)]./1000)
    
    %I - Ortho
    figure('Renderer', 'Painters')
    hold on
    plot(trn/1000,Apt_O(1,:),'LineWidth',2,'Color',[160/255 160/255 160/255])
    plot(trn/1000,Apt_O(3,:),'LineWidth',2,'Color',[10/255 206/255 10/255])
    plot(trn/1000,Apt_O(2,:),'LineWidth',2,'Color',[144/255 0 253/255])
    legend('No Drug','Bicuculine','Muscimol')
    yEnd=ylim;
    ylabel('PGC Firing Rate (Hz)')
    xlabel('Time (s)')
    title('Ortho')
    set(gca,'FontSize',18)
    set(gca,'XLim',[trn(1) trn(end)]./1000)
    %PGC - Retro
    figure('Renderer', 'Painters')
    hold on
    plot(trn/1000,Apt_R(1,:),'LineWidth',2,'Color',[160/255 160/255 160/255])
    plot(trn/1000,Apt_R(3,:),'LineWidth',2,'Color',[10/255 206/255 10/255])
    plot(trn/1000,Apt_R(2,:),'LineWidth',2,'Color',[144/255 0 253/255])
    legend('No Drug','Bicuculine','Muscimol')
    ylim(yEnd)
    ylabel('I Firing Rate (Hz)')
    xlabel('Time (s)')
    title('Retro')
    set(gca,'FontSize',18)
    set(gca,'XLim',[trn(1) trn(end)]./1000)
    
    %Plot weights and thresholds over time
    %Ortho
    figure('Renderer', 'Painters')
    hold on
    plot(trn/1000,w_MPt_O,'LineWidth',2)
    legend('wMP ND','wMP Bic','wMP Mus')
    xlabel('Time (s)')
    ylabel('Syn Weights','FontSize',12)
    title('Ortho','FontSize',18)
    
    %Retro
    figure('Renderer', 'Painters')
    hold on
    plot(trn/1000,w_MPt_R,'LineWidth',2)
    legend('wMP ND','wMP Bic','wMP Mus')
    xlabel('Time (s)')
    ylabel('Syn Weights','FontSize',12)
    title('Retro','FontSize',18)
    set(gca,'XLim',[trn(1) trn(end)]./1000)
end