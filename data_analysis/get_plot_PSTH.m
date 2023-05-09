% Script to get & plot PSTHs of RAW firing rate (within modality, across
% modality), heatmap of trial-avg firing rate and trial-std of firing,
% plots relating time-avg firing rate vs std dev (across trials), tables of correlation

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
        sRetTmp=sRET(FirstSpon:LastEvok,:,:);
        sOrTmp=sOR(FirstSpon:LastEvok,:,:);
        AllRat_sRET=cat(3,AllRat_sRET,sRetTmp);
        AllRat_sOR=cat(3,AllRat_sOR,sOrTmp);
    end %all individual rats
    DrugStruct{drug_to_keep+1,1}=AllRat_sOR;
    DrugStruct{drug_to_keep+1,2}=AllRat_sRET;
end %all drugs
PSTH_Orth_ND=DrugStruct{1,1}; %counts, NOT PSTH (divide by time below)
PSTH_Ret_ND=DrugStruct{1,2};
PSTH_Orth_Bic=DrugStruct{2,1};
PSTH_Ret_Bic=DrugStruct{2,2};
PSTH_Orth_Mus=DrugStruct{3,1};
PSTH_Ret_Mus=DrugStruct{3,2};
% calc PSTH, avg over trials
mnPSTH_Ret_ND=squeeze(mean(PSTH_Ret_ND,2)); mnPSTH_Or_ND=squeeze(mean(PSTH_Orth_ND,2));
mnPSTH_Ret_Bic=squeeze(mean(PSTH_Ret_Bic,2)); mnPSTH_Or_Bic=squeeze(mean(PSTH_Orth_Bic,2));
mnPSTH_Ret_Mus=squeeze(mean(PSTH_Ret_Mus,2)); mnPSTH_Or_Mus=squeeze(mean(PSTH_Orth_Mus,2));
% calc STD over trials
stPSTH_Ret_ND=squeeze(std(PSTH_Ret_ND,0,2)); stPSTH_Or_ND=squeeze(std(PSTH_Orth_ND,0,2));
stPSTH_Ret_Bic=squeeze(std(PSTH_Ret_Bic,0,2)); stPSTH_Or_Bic=squeeze(std(PSTH_Orth_Bic,0,2));
stPSTH_Ret_Mus=squeeze(std(PSTH_Ret_Mus,0,2)); stPSTH_Or_Mus=squeeze(std(PSTH_Orth_Mus,0,2));

%Plot
tme=(-2:0.1:2)'; %must match with lines 32,34; length is =LastEvok-FirstSpon+1

figure
subplot(3,1,1)
hold on
plot(tme,mean(mnPSTH_Or_ND,2)./setParam.Twin,'b','LineWidth',2)
plot(tme,mean(mnPSTH_Ret_ND,2)./setParam.Twin,'r','LineWidth',2)
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
title('ND')
set(gca,'FontSize',18)
axis([-.5 1 0 22])
hold off

subplot(3,1,2)
hold on
plot(tme,mean(mnPSTH_Or_Bic,2)./setParam.Twin,'b','LineWidth',2)
plot(tme,mean(mnPSTH_Ret_Bic,2)./setParam.Twin,'r','LineWidth',2)
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
title('Bic')
set(gca,'FontSize',18)
axis([-.5 1 0 22])
hold off

subplot(3,1,3)
hold on
plot(tme,mean(mnPSTH_Or_Mus,2)./setParam.Twin,'b','LineWidth',2)
plot(tme,mean(mnPSTH_Ret_Mus,2)./setParam.Twin,'r','LineWidth',2)
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
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
xlabel('Time (s)')
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
xlabel('Time (s)')
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
sr_ND_or=stPSTH_Or_ND./setParam.Twin; %making into rate (Hz)
sr_Bic_or=stPSTH_Or_Bic./setParam.Twin;
sr_Mus_or=stPSTH_Or_Mus./setParam.Twin;
sr_ND_ret=stPSTH_Ret_ND./setParam.Twin;
sr_Bic_ret=stPSTH_Ret_Bic./setParam.Twin;
sr_Mus_ret=stPSTH_Ret_Mus./setParam.Twin;
% % get time-avg mean rates; next set of Vars for DecAcc_mnGreater_raw.mat
% mean_nd_or=mean(fr_ND_or)';
% mean_nd_ret=mean(fr_ND_ret)';
% mean_bic_ret=mean(fr_Bic_ret)';
% mean_bic_or=mean(fr_Bic_or)';
% mean_mus_or=mean(fr_Mus_or)';
% mean_mus_ret=mean(fr_Mus_ret)';
% id_eq_Bic=(mean_bic_ret==mean_bic_or);
% id_eq_Mus=(mean_mus_ret==mean_mus_or);
% id_eq_ND=(mean_nd_ret==mean_nd_or);
% id_oGr_Bic=(mean_bic_or>mean_bic_ret & id_eq_Bic==0);
% id_oGr_Mus=(mean_mus_or>mean_mus_ret & id_eq_Mus==0);
% id_oGr_ND=(mean_nd_or>mean_nd_ret & id_eq_ND==0);
% id_rGr_Bic=(mean_bic_or<mean_bic_ret & id_eq_Bic==0);
% id_rGr_Mus=(mean_mus_or<mean_mus_ret & id_eq_Mus==0);
% id_rGr_ND=(mean_nd_or<mean_nd_ret & id_eq_ND==0);
% De_oGr_Bic=DecA_Bic(id_oGr_Bic);
% De_oGr_Mus=DecA_Mus(id_oGr_Mus);
% De_oGr_ND=DecA_ND(id_oGr_ND);
% De_rGo_Bic=DecA_Bic(id_rGr_Bic);
% De_rGo_Mus=DecA_Mus(id_rGr_Mus);
% De_rGo_ND=DecA_ND(id_rGr_Mus);

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

%% colormap plots whole population, mean
nmC_nd=size(fr_ND_or,2); nmC_bic=size(fr_Bic_or,2); nmC_mus=size(fr_Mus_or,2);
mnFR=min([fr_ND_or(:);fr_ND_ret(:);fr_Bic_or(:);fr_Bic_ret(:);fr_Mus_or(:);fr_Mus_ret(:)]);
mxFR=100;%max([fr_ND_or(:);fr_ND_ret(:);fr_Bic_or(:);fr_Bic_ret(:);fr_Mus_or(:);fr_Mus_ret(:)]);
% get scaled colormap
x = 1:512;
x = x-(0.9*mxFR-mnFR)*512/(mxFR-mnFR); %centerPoint=0.9*mxFR
x = 5*x/max(abs(x));  %scaleIntensity=5
x = sign(x).* exp(abs(x));
x = x - min(x); x = x*511/max(x)+1;
newMap = interp1(x, parula(512), 1:512);

[tmp,id1]=sort(mean(fr_ND_or));
figure('Renderer', 'painters', 'Position', [500 1500 650 300])
pcolor(tme,1:nmC_nd,fr_ND_or(:,id1)')
shading flat
colormap(newMap)
clim([mnFR mxFR])
%set(gca,'ColorScale','log')
colorbar
colormap parula
title('Mean firing rate, Ort-ND')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Cell Idnex (sorted)')
set(gca,'XLim',[-1 2])

[tmp,id2]=sort(mean(fr_ND_ret));
figure('Renderer', 'painters', 'Position', [500 1500 650 300])
pcolor(tme,1:nmC_nd,fr_ND_ret(:,id2)')
shading flat
colormap(newMap)
clim([mnFR mxFR])
colorbar
colormap parula
title('Mean firing rate, Ret-ND')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Cell Idnex (sorted)')
set(gca,'XLim',[-1 2])

[tmp,id3]=sort(mean(fr_Bic_or));
figure('Renderer', 'painters', 'Position', [500 1500 650 300])
pcolor(tme,1:nmC_bic,fr_Bic_or(:,id3)')
shading flat
colormap(newMap)
clim([mnFR mxFR])
colorbar
colormap parula
title('Mean firing rate, Ort-Bic')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Cell Idnex (sorted)')
set(gca,'XLim',[-1 2])

[tmp,id4]=sort(mean(fr_Bic_ret));
figure('Renderer', 'painters', 'Position', [500 1500 650 300])
pcolor(tme,1:nmC_bic,fr_Bic_ret(:,id4)')
shading flat
colormap(newMap)
clim([mnFR mxFR])
colorbar
colormap parula
title('Mean firing rate, Ret-Bic')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Cell Idnex (sorted)')
set(gca,'XLim',[-1 2])

[tmp,id5]=sort(mean(fr_Mus_or));
figure('Renderer', 'painters', 'Position', [500 1500 650 300])
pcolor(tme,1:nmC_mus,fr_Mus_or(:,id5)')
shading flat
colormap(newMap)
clim([mnFR mxFR])
colorbar
colormap parula
title('Mean firing rate, Ort-Mus')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Cell Idnex (sorted)')
set(gca,'XLim',[-1 2])

[tmp,id6]=sort(mean(fr_Mus_ret));
figure('Renderer', 'painters', 'Position', [500 1500 650 300])
pcolor(tme,1:nmC_mus,fr_Mus_ret(:,id6)')
shading flat
colormap(newMap)
clim([mnFR mxFR])
colorbar
colormap parula
title('Mean firing rate, Ret-Mus')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Cell Idnex (sorted)')
set(gca,'XLim',[-1 2])

%% colormap plots whole population, std
mnSR=min([sr_ND_or(:);sr_ND_ret(:);sr_Bic_or(:);sr_Bic_ret(:);sr_Mus_or(:);sr_Mus_ret(:)]);
mxSR=44.3844;%max([sr_ND_or(:);sr_ND_ret(:);sr_Bic_or(:);sr_Bic_ret(:);sr_Mus_or(:);sr_Mus_ret(:)]);
% get scaled colormap
x = 1:512;
x = x-(0.9*mxFR-mnFR)*512/(mxFR-mnFR); %centerPoint=0.9*mxFR
x = 1*x/max(abs(x));  %scaleIntensity=5
x = sign(x).* exp(abs(x));
x = x - min(x); x = x*511/max(x)+1;
newMapS = interp1(x, parula(512), 1:512);
%[tmp,id]=sort(mean(sr_ND_or));
figure('Renderer', 'painters', 'Position', [500 1500 650 300])
pcolor(tme,1:nmC_nd,sr_ND_or(:,id1)')
shading flat
colormap(newMapS)
clim([mnSR mxSR])
colorbar
colormap parula
title('STD firing rate, Ort-ND')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Cell Idnex (sorted)')
set(gca,'XLim',[-1 2])

%[tmp,id]=sort(mean(sr_ND_ret));
figure('Renderer', 'painters', 'Position', [500 1500 650 300])
pcolor(tme,1:nmC_nd,sr_ND_ret(:,id2)')
shading flat
colormap(newMapS)
clim([mnSR mxSR])
colorbar
colormap parula
title('STD firing rate, Ret-ND')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Cell Idnex (sorted)')
set(gca,'XLim',[-1 2])

%[tmp,id]=sort(mean(sr_Bic_or));
figure('Renderer', 'painters', 'Position', [500 1500 650 300])
pcolor(tme,1:nmC_bic,sr_Bic_or(:,id3)')
shading flat
colormap(newMapS)
clim([mnSR mxSR])
colorbar
colormap parula
title('STD firing rate, Ort-Bic')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Cell Idnex (sorted)')
set(gca,'XLim',[-1 2])

%[tmp,id]=sort(mean(sr_Bic_ret));
figure('Renderer', 'painters', 'Position', [500 1500 650 300])
pcolor(tme,1:nmC_bic,sr_Bic_ret(:,id4)')
shading flat
colormap(newMapS)
clim([mnSR mxSR])
colorbar
colormap parula
title('STD firing rate, Ret-Bic')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Cell Idnex (sorted)')
set(gca,'XLim',[-1 2])

%[tmp,id]=sort(mean(sr_Mus_or));
figure('Renderer', 'painters', 'Position', [500 1500 650 300])
pcolor(tme,1:nmC_mus,sr_Mus_or(:,id5)')
shading flat
colormap(newMapS)
clim([mnSR mxSR])
colorbar
colormap parula
title('STD firing rate, Ort-Mus')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Cell Idnex (sorted)')
set(gca,'XLim',[-1 2])

%[tmp,id]=sort(mean(sr_Mus_ret));
figure('Renderer', 'painters', 'Position', [500 1500 650 300])
pcolor(tme,1:nmC_mus,sr_Mus_ret(:,id6)')
shading flat
colormap(newMapS)
clim([mnSR mxSR])
colorbar
colormap parula
title('STD firing rate, Ret-Mus')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Cell Idnex (sorted)')
set(gca,'XLim',[-1 2])

%% plot time-avg evoked mean vs std for each MC/modality/regime, to see if there is a relationship
mn_ND_or=mean(fr_ND_or(21:end,:))'; %mean b/c fr_ & sr_ already in Hz
mn_ND_rt=mean(fr_ND_ret(21:end,:))';
mn_Bic_or=mean(fr_Bic_or(21:end,:))';
mn_Bic_rt=mean(fr_Bic_ret(21:end,:))';
mn_Mus_or=mean(fr_Mus_or(21:end,:))';
mn_Mus_rt=mean(fr_Mus_ret(21:end,:))';
std_ND_or=mean(sr_ND_or(21:end,:))';
std_ND_rt=mean(sr_ND_ret(21:end,:))';
std_Bic_or=mean(sr_Bic_or(21:end,:))';
std_Bic_rt=mean(sr_Bic_ret(21:end,:))';
std_Mus_or=mean(sr_Mus_or(21:end,:))';
std_Mus_rt=mean(sr_Mus_ret(21:end,:))';
figure
hold on
plot(mn_ND_or,std_ND_or,'.','MarkerSize',18,'color',.502*ones(1,3))
plot(mn_Bic_or,std_Bic_or,'.','MarkerSize',18,'color',[0 200 0]./255)
plot(mn_Mus_or,std_Mus_or,'.','MarkerSize',18,'color',[140 0 255]./255)
set(gca,'FontSize',20)
xlabel('Time-Averaged Evoked Firing Rate (Hz)')
ylabel('Time-Averaged STD of Evoked Firing Rate (Hz)')
axis([0 150 0 60])
figure
hold on
plot(mn_ND_rt,std_ND_rt,'.','MarkerSize',18,'color',.502*ones(1,3))
plot(mn_Bic_rt,std_Bic_rt,'.','MarkerSize',18,'color',[0 200 0]./255)
plot(mn_Mus_rt,std_Mus_rt,'.','MarkerSize',18,'color',[140 0 255]./255)
set(gca,'FontSize',20)
xlabel('Time-Averaged Evoked Firing Rate (Hz)')
ylabel('Time-Averaged STD of Evoked Firing Rate (Hz)')
axis([0 150 0 60])

drgPrep={'ND';'Bic';'Mus'}; 

%showing correlations betwn mean & std
pearsonsCorr=[corr(mn_ND_or,std_ND_or) ; corr(mn_Bic_or,std_Bic_or) ; corr(mn_Mus_or,std_Mus_or)];
RankCorr=[corr(mn_ND_or,std_ND_or,'type','Spearman') ; corr(mn_Bic_or,std_Bic_or,'type','Spearman')  ; corr(mn_Mus_or,std_Mus_or,'type','Spearman') ];
T_ort=table(drgPrep,pearsonsCorr,RankCorr)

pearsonsCorr=[corr(mn_ND_rt,std_ND_rt) ; corr(mn_Bic_rt,std_Bic_rt) ; corr(mn_Mus_rt,std_Mus_rt)];
RankCorr=[corr(mn_ND_rt,std_ND_rt,'type','Spearman') ; corr(mn_Bic_rt,std_Bic_rt,'type','Spearman')  ; corr(mn_Mus_rt,std_Mus_rt,'type','Spearman') ];
T_ret=table(drgPrep,pearsonsCorr,RankCorr)


%% corresponding graph of indiv PSTH
if(0) %same info as heatmaps, not useful
figure('Renderer', 'Painters');
for j=nmC_nd:-1:1
    jFT=ceil(size(newMap,1)*j/nmC_nd);
    semilogy(tme,fr_ND_or(:,id1(j)),'color',newMap(jFT,:))
    if(j==nmC_nd) 
        hold on
    end
    alpha((nmC_nd-j)/(nmC_nd-1))
end
title('Mean firing rate, Ort-ND')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
axis([-1 2 0 250])
figure('Renderer', 'Painters');
for j=nmC_nd:-1:1
    jFT=ceil(size(newMap,1)*j/nmC_nd);
    semilogy(tme,fr_ND_ret(:,id2(j)),'color',newMap(jFT,:))
    if(j==nmC_nd)
        hold on
    end
    alpha((nmC_nd-j)/(nmC_nd-1))
end
title('Mean firing rate, Retr-ND')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
axis([-1 2 0 250])

figure('Renderer', 'Painters');
for j=nmC_bic:-1:1
    jFT=ceil(size(newMap,1)*j/nmC_bic);
    semilogy(tme,fr_Bic_or(:,id3(j)),'color',newMap(jFT,:))
    if(j==nmC_bic)
        hold on
    end
    alpha((nmC_bic-j)/(nmC_bic-1))
end
title('Mean firing rate, Ort-Bic')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
axis([-1 2 0 250])
figure('Renderer', 'Painters');
for j=nmC_bic:-1:1
    jFT=ceil(size(newMap,1)*j/nmC_bic);
    semilogy(tme,fr_Bic_ret(:,id4(j)),'color',newMap(jFT,:))
    if(j==nmC_bic)
        hold on
    end
    alpha((nmC_bic-j)/(nmC_bic-1))
end
title('Mean firing rate, Retr-Bic')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
axis([-1 2 0 250])

figure('Renderer', 'Painters');
for j=nmC_mus:-1:1
    jFT=ceil(size(newMap,1)*j/nmC_mus);
    semilogy(tme,fr_Mus_or(:,id5(j)),'color',newMap(jFT,:))
    if(j==nmC_mus)
        hold on
    end
    alpha((nmC_mus-j)/(nmC_mus-1))
end
title('Mean firing rate, Ort-Mus')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
axis([-1 2 0 250])
figure('Renderer', 'Painters');
for j=nmC_mus:-1:1
    jFT=ceil(size(newMap,1)*j/nmC_mus);
    semilogy(tme,fr_Mus_ret(:,id6(j)),'color',newMap(jFT,:))
    if(j==nmC_mus)
        hold on
    end
    alpha((nmC_mus-j)/(nmC_mus-1))
end
title('Mean firing rate, Retr-Mus')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
axis([-1 2 0 250])

% --- repeat for STD----
figure('Renderer', 'Painters');
for j=nmC_nd:-1:1
    jFT=ceil(size(newMap,1)*j/nmC_nd);
    semilogy(tme,sr_ND_or(:,id1(j)),'color',newMapS(jFT,:))
    if(j==nmC_nd)
        hold on
    end
    alpha((nmC_nd-j)/(nmC_nd-1))
end
title('STD firing rate, Ort-ND')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('STD FR (Hz)')
axis([-1 2 0 105])
figure('Renderer', 'Painters');
for j=nmC_nd:-1:1
    jFT=ceil(size(newMap,1)*j/nmC_nd);
    semilogy(tme,sr_ND_ret(:,id2(j)),'color',newMapS(jFT,:))
    if(j==nmC_nd)
        hold on
    end
    alpha((nmC_nd-j)/(nmC_nd-1))
end
title('STD firing rate, Retr-ND')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('STD FR (Hz)')
axis([-1 2 0 105])

figure('Renderer', 'Painters');
for j=nmC_bic:-1:1
    jFT=ceil(size(newMap,1)*j/nmC_bic);
    semilogy(tme,sr_Bic_or(:,id3(j)),'color',newMapS(jFT,:))
    if(j==nmC_bic)
        hold on
    end
    alpha((nmC_bic-j)/(nmC_bic-1))
end
title('STD FR, Ort-Bic')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('STD FR (Hz)')
axis([-1 2 0 105])
figure('Renderer', 'Painters')
for j=nmC_bic:-1:1
    jFT=ceil(size(newMap,1)*j/nmC_bic);
    semilogy(tme,sr_Bic_ret(:,id4(j)),'color',newMapS(jFT,:))
    if(j==nmC_bic)
        hold on
    end
    alpha((nmC_bic-j)/(nmC_bic-1))
end
title('STD FR, Retr-Bic')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('STD FR (Hz)')
axis([-1 2 0 105])

figure('Renderer', 'Painters');
for j=nmC_mus:-1:1
    jFT=ceil(size(newMap,1)*j/nmC_mus);
    semilogy(tme,sr_Mus_or(:,id5(j)),'color',newMapS(jFT,:))
    if(j==nmC_mus)
        hold on
    end
    alpha((nmC_mus-j)/(nmC_mus-1))
end
title('STD FR, Ort-Mus')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('STD FR (Hz)')
axis([-1 2 0 105])
figure('Renderer', 'Painters');
for j=nmC_mus:-1:1
    jFT=ceil(size(newMap,1)*j/nmC_mus);
    semilogy(tme,sr_Mus_ret(:,id6(j)),'color',newMapS(jFT,:))
    if(j==nmC_mus)
        hold on
    end
    alpha((nmC_mus-j)/(nmC_mus-1))
end
title('STD FR, Retr-Mus')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('STD FR (Hz)')
axis([-1 2 0 105])
end