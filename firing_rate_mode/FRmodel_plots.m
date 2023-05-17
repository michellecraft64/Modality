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
TimeVars=struct('dt',0.01,'lenSpon',500,'lenEvok',2000,'lenTrn',1000); %opt Twin = 900ms evoked (2s spont)

% % weights below & lamO/lamR 
%with original lamO/lamR
alphMax_or_MC=80; %Constant multiplier for input to MC
alphMax_or_I=95; %Constant multiplier for input to I
alphMax_ret_MC=50; %Constant multiplier for input to MC
alphMax_ret_I=15; %Constant multiplier for input to I

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
lamO=struct('tsft',0,'lmOsp',.14,'lmOevk',0.1035,'tauO',95,'tauO2',100,'ssV',.15); % based on Craft et al 21
%tmShift=1; %in ms
TimeVars.lenEvokS=tmeEvok;%[tmeEvok; tmeEvok(end)+(dt:dt:tmShift)']; %take nuevok out further by tmShift
%nmShift=length((dt:dt:tmShift)); %# elem must shift by
nuevokO=lamO.lmOevk*(-(TimeVars.lenEvokS+lamO.tsft).*exp(-(TimeVars.lenEvokS+lamO.tsft)/lamO.tauO) + (TimeVars.lenEvokS+lamO.tsft).*exp(-(TimeVars.lenEvokS+lamO.tsft)/lamO.tauO2))+lamO.ssV;
%nuevokO=nuevokO(nmShift+1:end); %shift over by tmShift
lam_O=[lamO.lmOsp*ones(length(tmeSpon)-1,1);nuevokO];
%Retro input
lamR=struct('tsft',200,'lmRspon',.14,'lmRevk',0.001,'tau_R',5,'tauR2',500,'ssV',.25); % based on Craft et al 21
nuevok_R=lamR.lmRevk*(-(TimeVars.lenEvokS+lamR.tsft).*exp(-(TimeVars.lenEvokS+lamR.tsft)/lamR.tau_R) + (TimeVars.lenEvokS+lamR.tsft).*exp(-(TimeVars.lenEvokS+lamR.tsft)/lamR.tauR2))+lamR.ssV;
idSplt=round(202.2/dt+1); %see where plot(Tmev,lam_O,Tmev,lam_R,'r') intersects
nuevok_R(1:idSplt)=nuevokO(1:idSplt); %split, first part of Retro like Ortho
lam_R=[lamR.lmRspon*ones(length(tmeSpon)-1,1);nuevok_R];
%plot(Tmev,lam_O,Tmev,lam_R,'r')

InputCurrent_Or=struct('Mc',lam_O*alphMax_or_MC,'Pg',lam_O*alphMax_or_I);
InputCurrent_Ret=struct('Mc',lam_R*alphMax_ret_MC,'Pg',lam_R*alphMax_ret_I);

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

%Iterate through params
for i=1:numParams
    %Define specific time constants and coupling strengths
    TimeConstants=struct('tau_M',10,'tau_P',5.5,'tdE',ParamMat(1,i),'trE',ParamMat(2,i),'tdI',10,'trI',2,'tw',tw);
    % td[E/I]=decay time of syn, tr[E/I]=rise time of syn
    CpldWx=struct('w_I',1,'w_P',ParamMat(3,i),'wmp0',1); 
    
    switch flag_Regime
        case 1
            %call for Ortho
            [Am_O(i,:),Ap_O(i,:),syEo(i,:),syIo(i,:),w_MP_O(i,:),Amt_O(i,:),Apt_O(i,:),w_MPt_O(i,:)]=WC_EIcells(CpldWx,TimeVars,TimeConstants,xMin,InputCurrent_Or,crr);
            %call for Ret
            [Am_R(i,:),Ap_R(i,:),syEr(i,:),syIr(i,:),w_MP_R(i,:),Amt_R(i,:),Apt_R(i,:),w_MPt_R(i,:)]=WC_EIcells(CpldWx,TimeVars,TimeConstants,xMin,InputCurrent_Ret,crr);
        case 2
            %call for Ortho
            [Am_O(i,:),Ap_O(i,:),syEo(i,:),syIo(i,:),w_MP_O(i,:),Amt_O(i,:),Apt_O(i,:),w_MPt_O(i,:)]=WC_EIcellsF(CpldWx,TimeVars,TimeConstants,xMin,InputCurrent_Or,crr);
            %call for Ret
            [Am_R(i,:),Ap_R(i,:),syEr(i,:),syIr(i,:),w_MP_R(i,:),Amt_R(i,:),Apt_R(i,:),w_MPt_R(i,:)]=WC_EIcellsF(CpldWx,TimeVars,TimeConstants,xMin,InputCurrent_Ret,crr);
    end

end

% Am_base_O=Am_O;
% Am_base_R=Am_R;
% Ap_base_O=Ap_O;
% Ap_base_R=Ap_R;
% save dBase_noDepr Tmev Am_base_O Am_base_R Ap_base_O Ap_base_R alphMax*

%% Neg bin decode acc
aspPv=[0.32 0.32; .32 .052; .32 .53];
idZer=round((0-Tmev(1))/dt*1000+1); idEvD=round((0.9-Tmev(1))/dt*1000+1);
if(NBdAcc_Flag==1)
    tic
    Ntr=50000; %aspirational p for neg binomial, p=1 is Poisson, p inverse prop to variance
    AccMat=zeros(numParams,4);
    for k=1:numParams
        
        nuO=squeeze(Am_O(k,idZer:idEvD));
        % Negative Binomial
        meanOc=sum(nuO)*dt/(1000*0.9); %b/c Tmev in s & nu0 is per ms, convt to Hz!
        sc(k,1)=aspPv(k,1)/(1-aspPv(k,1)); %r=sc*mu, sc=p/(1-p) in neg binomials
        rO=meanOc*sc(k,1);
        
        nuR=squeeze(Am_R(k,idZer:idEvD));
        meanRc=sum(nuR)*dt/(1000*0.9); %b/c Tmev in s & nu0 is per msm convt to Hz!
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

    figure
    hold on
    plot(Tmev,InputCurrent_Or.Mc,'b','LineWidth',2)
    plot(Tmev,InputCurrent_Ret.Mc,'r','LineWidth',2)
    plot(Tmev,InputCurrent_Or.Pg,'c-','LineWidth',1)
    plot(Tmev,InputCurrent_Ret.Pg,'m-','LineWidth',1)
    set(gca,'FontSize',18)
    xlabel('Time (s)')
    ylabel('ORN input (a.u.)')

%Plot different "drugs"
load dExperim_FR.mat
%MC - Ortho
ccD=[128.01 128.01 128.01; 140 0 255; 0 200 0]./255;
figure('Renderer', 'Painters')
hold on
plot(Tmev,Am_O(1,:),'LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Am_O(3,:),'LineWidth',2,'Color',ccD(3,:))
plot(Tmev,Am_O(2,:),'LineWidth',2,'Color',ccD(2,:))
%legend('No Drug','Bicuculine','Muscimol')
yEnd=ylim;
ylabel('MC Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'FontSize',18)
set(gca,'XLim',[Tmev(1) Tmev(end)])
load dBase_noDepr.mat
plot(Tmev,Am_base_O(1,:),'--','LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Am_base_O(2,:),'--','LineWidth',2,'Color',ccD(2,:))
plot(Tmev,Am_base_O(3,:),'--','LineWidth',2,'Color',ccD(3,:))
%plot w/ Experim Data:
plot(tme,paFR_PSTH_Or_ND,'.-','LineWidth',.25,'Color',[160/255 160/255 160/255],'MarkerSize',22)
plot(tme,paFR_PSTH_Or_Bic,'.-','LineWidth',.25,'Color',[10/255 206/255 10/255],'MarkerSize',22)
plot(tme,paFR_PSTH_Or_Mus,'.-','LineWidth',.25,'Color',[144/255 0 253/255],'MarkerSize',22)
axis([-.5 2 0 25])

%MC - Retro
figure('Renderer', 'Painters')
hold on
plot(Tmev,Am_R(1,:),'LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Am_R(3,:),'LineWidth',2,'Color',ccD(3,:))
plot(Tmev,Am_R(2,:),'LineWidth',2,'Color',ccD(2,:))
%legend('No Drug','Bicuculine','Muscimol')
ylim(yEnd)
ylabel('MC Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'FontSize',18)
set(gca,'XLim',[Tmev(1) Tmev(end)])
plot(Tmev,Am_base_R(1,:),'--','LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Am_base_R(2,:),'--','LineWidth',2,'Color',ccD(2,:))
plot(Tmev,Am_base_R(3,:),'--','LineWidth',2,'Color',ccD(3,:))
% plot with Experim Data:
plot(tme,paFR_PSTH_Ret_ND,'.-','LineWidth',.25,'Color',[160/255 160/255 160/255],'MarkerSize',22)
plot(tme,paFR_PSTH_Ret_Bic,'.-','LineWidth',.25,'Color',[10/255 206/255 10/255],'MarkerSize',22)
plot(tme,paFR_PSTH_Ret_Mus,'.-','LineWidth',.25,'Color',[144/255 0 253/255],'MarkerSize',22)
axis([-.5 2 0 25])

if(0) %for the refs
    figure
    subplot(3,1,1) 
    hold on
    plot(Tmev,Am_O(1,:),'b','LineWidth',2)
    plot(Tmev,Am_R(1,:),'r','LineWidth',2)
    plot(tme,paFR_PSTH_Or_ND,'b.-','LineWidth',1,'MarkerSize',12)
    plot(tme,paFR_PSTH_Ret_ND,'r.-','LineWidth',1,'MarkerSize',12)
    axis([-.5 2 0 25])
    subplot(3,1,2) 
    hold on
    plot(Tmev,Am_O(3,:),'b','LineWidth',2)
    plot(Tmev,Am_R(3,:),'r','LineWidth',2)
    plot(tme,paFR_PSTH_Or_Bic,'b.-','LineWidth',1,'MarkerSize',12)
    plot(tme,paFR_PSTH_Ret_Bic,'r.-','LineWidth',1,'MarkerSize',12)
    axis([-.5 2 0 25])
    subplot(3,1,3) 
    hold on
    plot(Tmev,Am_O(2,:),'b','LineWidth',2)
    plot(Tmev,Am_R(2,:),'r','LineWidth',2)
    plot(tme,paFR_PSTH_Or_Mus,'b.-','LineWidth',1,'MarkerSize',12)
    plot(tme,paFR_PSTH_Ret_Mus,'r.-','LineWidth',1,'MarkerSize',12)
    axis([-.5 2 0 25])
end

%I - Ortho
figure('Renderer', 'Painters')
hold on
plot(Tmev,Ap_O(1,:),'LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Ap_O(3,:),'LineWidth',2,'Color',ccD(3,:))
plot(Tmev,Ap_O(2,:),'LineWidth',2,'Color',ccD(2,:))
legend('No Drug','Bicuculine','Muscimol')
yEnd=ylim;
ylabel('Ortho I Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'FontSize',18)
set(gca,'XLim',[Tmev(1) Tmev(end)])
plot(Tmev,Ap_base_O(1,:),'--','LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Ap_base_O(2,:),'--','LineWidth',2,'Color',ccD(2,:))
plot(Tmev,Ap_base_O(3,:),'--','LineWidth',2,'Color',ccD(3,:))
axis([-.5 2 0 70])

%I - Retro
figure('Renderer', 'Painters')
hold on
plot(Tmev,Ap_R(1,:),'LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Ap_R(3,:),'LineWidth',2,'Color',ccD(3,:))
plot(Tmev,Ap_R(2,:),'LineWidth',2,'Color',ccD(2,:))
legend('No Drug','Bicuculine','Muscimol')
ylim(yEnd)
ylabel('Retro I Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'FontSize',18)
set(gca,'XLim',[Tmev(1) Tmev(end)])
plot(Tmev,Ap_base_R(1,:),'--','LineWidth',2,'Color',ccD(1,:))
plot(Tmev,Ap_base_R(2,:),'--','LineWidth',2,'Color',ccD(2,:))
plot(Tmev,Ap_base_R(3,:),'--','LineWidth',2,'Color',ccD(3,:))
axis([-.5 2 0 70])

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

%% table of relative change
if(1)
    idZer=round((0-Tmev(1))/dt*1000+1); idEvD=round((0.9-Tmev(1))/dt*1000+1); 
or_nd_Dep=sum(Am_O(1,idZer:idEvD))*dt/(1000*0.9);
or_bic_Dep=sum(Am_O(3,idZer:idEvD))*dt/(1000*0.9);
or_mus_Dep=sum(Am_O(2,idZer:idEvD))*dt/(1000*0.9);
    rt_nd_Dep=sum(Am_R(1,idZer:idEvD))*dt/(1000*0.9);
    rt_bic_Dep=sum(Am_R(3,idZer:idEvD))*dt/(1000*0.9);
    rt_mus_Dep=sum(Am_R(2,idZer:idEvD))*dt/(1000*0.9);
or_nd_base=sum(Am_base_O(1,idZer:idEvD))*dt/(1000*0.9);
or_bic_base=sum(Am_base_O(3,idZer:idEvD))*dt/(1000*0.9);
or_mus_base=sum(Am_base_O(2,idZer:idEvD))*dt/(1000*0.9);
    rt_nd_base=sum(Am_base_R(1,idZer:idEvD))*dt/(1000*0.9);
    rt_bic_base=sum(Am_base_R(3,idZer:idEvD))*dt/(1000*0.9);
    rt_mus_base=sum(Am_base_R(2,idZer:idEvD))*dt/(1000*0.9);

drgComp={'ND/Bic';'ND/Mus'}; 
orth_NoDeprs=[or_bic_base-or_nd_base ; or_nd_base-or_mus_base  ]./or_nd_base;
orth_Deprs=[or_bic_Dep-or_nd_Dep ; or_nd_Dep-or_mus_Dep  ]./or_nd_Dep;
ret_NoDeprs=[rt_bic_base-rt_nd_base ; rt_nd_base-rt_mus_base  ]./rt_nd_base;
ret_Deprs=[rt_bic_Dep-rt_nd_Dep ; rt_nd_Dep-rt_mus_Dep  ]./rt_nd_Dep;
T_relChange=table(drgComp,orth_NoDeprs,orth_Deprs,ret_NoDeprs,ret_Deprs)

drgs={'ND';'Bic';'Mus'};
T_integFR=table(drgs,[or_nd_base; or_bic_base; or_mus_base],[or_nd_Dep; or_bic_Dep; or_mus_Dep],...
    [rt_nd_base; rt_bic_base; rt_mus_base],[rt_nd_Dep; rt_bic_Dep; rt_mus_Dep])

% dtD=tme(2)-tme(1);
% idZ_d=round((0-tme(1))/dtD+1); idEv_d=round((0.9-tme(1))/dtD+1); 
% or_nd_data=sum(paFR_PSTH_Or_ND(idZ_d:idEv_d))*dtD/(0.9);
% or_bic_data=sum(paFR_PSTH_Or_Bic(idZ_d:idEv_d))*dtD/(0.9);
% or_mus_data=sum(paFR_PSTH_Or_Mus(idZ_d:idEv_d))*dtD/(0.9);
%     rt_nd_data=sum(paFR_PSTH_Ret_ND(idZ_d:idEv_d))*dtD/(0.9);
%     rt_bic_data=sum(paFR_PSTH_Ret_Bic(idZ_d:idEv_d))*dtD/(0.9);
%     rt_mus_data=sum(paFR_PSTH_Ret_Mus(idZ_d:idEv_d))*dtD/(0.9);
% 
% ort_rel=[or_bic_data-or_nd_data ; or_nd_data-or_mus_data  ]./or_nd_data;
% ret_rel=[rt_bic_data-rt_nd_data ; rt_nd_data-rt_mus_data  ]./rt_nd_data;
% T_relChngData=table(drgComp,ort_rel,ret_rel)

%% table in manuscript for simulated trial-var (model) = mu/rho

drgComp={'Ortho';'Retro'}; 
Muscimol=[or_mus_Dep; rt_mus_Dep]./(aspPv(2,:)'); %recall 2nd row is Mus in this case
NoDrug=[or_nd_Dep; rt_nd_Dep]./(aspPv(1,:)');
Bicuculline=[or_bic_Dep; rt_bic_Dep]./(aspPv(3,:)');
T_trlVar=table(drgComp,Muscimol,NoDrug,Bicuculline)
end

end

