% Script to get & plot PSTHs
clear
% Parameters
odorName='EB';

% Initialize structure
DrugStruct=cell(3,2); %Rows 1) Or 2) Ret / Columns 1) ND 2) Bic 3) Mus
% Loop over all drugs
for drug_to_keep = 0:2 % DRUG: no drug = 0; Bicu = 1 (less GABA_a inhib); Musc = 2 (more GABA_a inhib)
    %Load good rat data
    if drug_to_keep==0
        ind=[1 2 6 7 8 9 10 11]; %actual good rats
        drugName='NoDrug';
    elseif drug_to_keep==1
        ind=[1 2 8 11]; %actual good rats
        drugName='Bic';
    elseif drug_to_keep == 2
        ind=[6 9 10]; %actual good rats
        drugName='Mus';
    end
    numRats=length(ind);
    AllRat_sRET=[];
    AllRat_sOR=[];
    for i=ind
        fileName=sprintf('Rat%d_IndCell_%s_%s.mat',i,odorName,drugName);
        %Load file, set parameters
        load(fileName) %Ortho/Retro raw spike counts (lenTime, numTrials, nID)
        numTrials=size(sOR,2)+size(sRET,2); %per or/ret cell
        numEvok=TimeVars.numEvok;
        StimShift=TimeVars.StimShift;
        FirstEvok=TimeVars.FirstEvok;
        LastEvok=FirstEvok+19; %TimeVars.LastEvok; %HARD CODED for opt twin-NEEDS CHANGING IF TWIN CHANGED
        LastSpon=TimeVars.LastSpon;
        FirstSpon=LastSpon-20; %TimeVars.FirstSpon; %HARD CODED for opt twin-NEEDS CHANGING IF TWIN CHANGED
        lenEvok=length(FirstEvok:1:LastEvok);
        lenSpon=length(FirstSpon:1:LastSpon);
        nOB=size(sOR,3);
        if nOB~=size(sRET,3)
            keyboard
        end
        sRetTmp=sRET(FirstSpon:LastEvok,:,:);
        sOrTmp=sOR(FirstSpon:LastEvok,:,:);
        AllRat_sRET=cat(3,AllRat_sRET,sRetTmp);
        AllRat_sOR=cat(3,AllRat_sOR,sOrTmp);
    end %all individual rats
    DrugStruct{drug_to_keep+1,1}=AllRat_sOR;
    DrugStruct{drug_to_keep+1,2}=AllRat_sRET;
end %all drugs
PSTH_Orth_ND=DrugStruct{1,1};
PSTH_Ret_ND=DrugStruct{1,2};
PSTH_Orth_Bic=DrugStruct{2,1};
PSTH_Ret_Bic=DrugStruct{2,2};
PSTH_Orth_Mus=DrugStruct{3,1};
PSTH_Ret_Mus=DrugStruct{3,2};
mnPSTH_Ret_ND=squeeze(mean(PSTH_Ret_ND,2)); mnPSTH_Or_ND=squeeze(mean(PSTH_Orth_ND,2));
mnPSTH_Ret_Bic=squeeze(mean(PSTH_Ret_Bic,2)); mnPSTH_Or_Bic=squeeze(mean(PSTH_Orth_Bic,2));
mnPSTH_Ret_Mus=squeeze(mean(PSTH_Ret_Mus,2)); mnPSTH_Or_Mus=squeeze(mean(PSTH_Orth_Mus,2));

%Plot
tme=(-2:0.1:2)'; %must match with lines 32,34; length is =LastEvok-FirstSpon+1

figure
subplot(3,1,1)
hold on
plot(tme,mean(mnPSTH_Or_ND,2)./setParam.Twin,'b','LineWidth',2)
plot(tme,mean(mnPSTH_Ret_ND,2)./setParam.Twin,'r','LineWidth',2)
ylabel('Firing Rate (Hz)')
xlabel('Time (ms)')
title('ND')
set(gca,'FontSize',18)
axis([-.5 1 0 22])
hold off

subplot(3,1,2)
hold on
plot(tme,mean(mnPSTH_Or_Bic,2)./setParam.Twin,'b','LineWidth',2)
plot(tme,mean(mnPSTH_Ret_Bic,2)./setParam.Twin,'r','LineWidth',2)
ylabel('Firing Rate (Hz)')
xlabel('Time (ms)')
title('Bic')
set(gca,'FontSize',18)
axis([-.5 1 0 22])
hold off

subplot(3,1,3)
hold on
plot(tme,mean(mnPSTH_Or_Mus,2)./setParam.Twin,'b','LineWidth',2)
plot(tme,mean(mnPSTH_Ret_Mus,2)./setParam.Twin,'r','LineWidth',2)
ylabel('Firing Rate (Hz)')
xlabel('Time (ms)')
title('Mus')
set(gca,'FontSize',18)
axis([-.5 1 0 22])
hold off

figure
hold on
plot(tme,mean(mnPSTH_Or_ND,2)./setParam.Twin,'LineWidth',2,'Color',[160/255 160/255 160/255])
plot(tme,mean(mnPSTH_Or_Bic,2)./setParam.Twin,'LineWidth',2,'Color',[10/255 206/255 10/255])
plot(tme,mean(mnPSTH_Or_Mus,2)./setParam.Twin,'LineWidth',2,'Color',[144/255 0 253/255])
legend('No Drug','Bicuculine','Muscimol','Location','northwest')
ylabel('Firing Rate (Hz)')
xlabel('Time (ms)')
title('Ortho')
set(gca,'FontSize',18)
axis([-2 2 0 22])
hold off

figure
hold on
plot(tme,mean(mnPSTH_Ret_ND,2)./setParam.Twin,'LineWidth',2,'Color',[160/255 160/255 160/255])
plot(tme,mean(mnPSTH_Ret_Bic,2)./setParam.Twin,'LineWidth',2,'Color',[10/255 206/255 10/255])
plot(tme,mean(mnPSTH_Ret_Mus,2)./setParam.Twin,'LineWidth',2,'Color',[144/255 0 253/255])
legend('No Drug','Bicuculine','Muscimol','Location','northwest')
ylabel('Firing Rate (Hz)')
xlabel('Time (ms)')
title('Retro')
set(gca,'FontSize',18)
axis([-2 2 0 22])
hold off

%% plotting p-values from t-test time series
alph=0.01; %signif value
hghtStr=0.2; %height of star in plot
LnmW=LastEvok-FirstSpon+1; %size of time (-2 to 2), etc
pVl_NDBic_or=zeros(LnmW,1);
pVl_NDMus_or=zeros(LnmW,1);
pVl_BicMus_or=zeros(LnmW,1);
pVl_NDBic_ret=zeros(LnmW,1);
pVl_NDMus_ret=zeros(LnmW,1);
pVl_BicMus_ret=zeros(LnmW,1);

fr_ND_or=mnPSTH_Or_ND./setParam.Twin;
fr_Bic_or=mnPSTH_Or_Bic./setParam.Twin;
fr_Mus_or=mnPSTH_Or_Mus./setParam.Twin;
fr_ND_ret=mnPSTH_Ret_ND./setParam.Twin;
fr_Bic_ret=mnPSTH_Ret_Bic./setParam.Twin;
fr_Mus_ret=mnPSTH_Ret_Mus./setParam.Twin;
for j=1:LnmW
    [~,p]=ttest2(fr_ND_or(j,:),fr_Bic_or(j,:),'VarType','unequal');
    pVl_NDBic_or(j)=p;
    [~,p]=ttest2(fr_ND_or(j,:),fr_Mus_or(j,:),'VarType','unequal');
    pVl_NDMus_or(j)=p;
    [~,p]=ttest2(fr_Bic_or(j,:),fr_Mus_or(j,:),'VarType','unequal');
    pVl_BicMus_or(j)=p;
    
    [~,p]=ttest2(fr_ND_ret(j,:),fr_Bic_ret(j,:),'VarType','unequal');
    pVl_NDBic_ret(j)=p;
    [~,p]=ttest2(fr_ND_ret(j,:),fr_Mus_ret(j,:),'VarType','unequal');
    pVl_NDMus_ret(j)=p;
    [~,p]=ttest2(fr_Bic_ret(j,:),fr_Mus_ret(j,:),'VarType','unequal');
    pVl_BicMus_ret(j)=p;
end

%for ortho
figure
subplot(1,3,1)
hold on
plot(tme,pVl_NDMus_or,'k')
plot(tme(pVl_NDMus_or<alph & tme>0),hghtStr*ones(size(tme(pVl_NDMus_or<alph & tme>0))),'k*')
plot(tme,alph*ones(size(tme)),'k--')
axis([-2 2 0 .21])
subplot(1,3,2)
hold on
plot(tme,pVl_NDBic_or,'k')
plot(tme(pVl_NDBic_or<alph & tme>0),hghtStr*ones(size(tme(pVl_NDBic_or<alph & tme>0))),'k*')
plot(tme,alph*ones(size(tme)),'k--')
axis([-2 2 0 .21])
subplot(1,3,3)
hold on
plot(tme,pVl_BicMus_or,'k')
plot(tme(pVl_BicMus_or<alph & tme>0),hghtStr*ones(size(tme(pVl_BicMus_or<alph & tme>0))),'k*')
plot(tme,alph*ones(size(tme)),'k--')
axis([-2 2 0 .21])

%for retro
figure
subplot(1,3,1)
hold on
plot(tme,pVl_NDMus_ret,'k')
plot(tme,alph*ones(size(tme)),'k--')
axis([-2 2 0 1])
subplot(1,3,2)
hold on
plot(tme,pVl_NDBic_ret,'k')
plot(tme,alph*ones(size(tme)),'k--')
axis([-2 2 0 1])
subplot(1,3,3)
hold on
plot(tme,pVl_BicMus_ret,'k')
plot(tme,alph*ones(size(tme)),'k--')
axis([-2 2 0 1])

