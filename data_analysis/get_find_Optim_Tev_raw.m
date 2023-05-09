%Script to plot and find Optimal Twin (t*), ONLY evoked b/c raw FR
clear
% Select parameters
odor_to_keep = 2; 
%For generating labels
or='Ortho';
ret='Retro';
odorName='EB';
sponTime=0;

fileName=sprintf('IndCellRAW_OptWin_%s_%dsSpon.mat',odorName,sponTime);

p_NDbic_0Spon=[];
p_NDmus_0Spon=[];
pW_NDbic_0Spon=[];
pW_NDmus_0Spon=[];

for j=1:10
    PlotTtl=sprintf('All rat net accuracies - %s %dms',odorName,j*100);
    NDacc=[];
    BICacc=[];
    MUSacc=[];
    load(fileName,'optThrsh_ND','optThrsh_Bic','optThrsh_Mus')
    NDnumRats=size(optThrsh_ND,1);
    BICnumRats=size(optThrsh_Bic,1);
    MUSnumRats=size(optThrsh_Mus,1);
    for i=1:NDnumRats
        NDall=optThrsh_ND{i,4};
        NDacc=[NDacc NDall(j,:)];
    end
    for i=1:BICnumRats
        BICall=optThrsh_Bic{i,4};
        BICacc=[BICacc BICall(j,:)];
    end
    for i=1:MUSnumRats
        MUSall=optThrsh_Mus{i,4};
        MUSacc=[MUSacc MUSall(j,:)];
    end
    [~,p_NDbic]=ttest2(NDacc,BICacc,'VarType','unequal');
    [~,p_NDmus]=ttest2(NDacc,MUSacc,'VarType','unequal');
    pW_NDbic=ranksum(NDacc,BICacc);
    pW_NDmus=ranksum(NDacc,MUSacc);
    
    p_NDbic_0Spon=[p_NDbic_0Spon p_NDbic];
    p_NDmus_0Spon=[p_NDmus_0Spon p_NDmus];
    pW_NDbic_0Spon=[pW_NDbic_0Spon pW_NDbic];
    pW_NDmus_0Spon=[pW_NDmus_0Spon pW_NDmus];

    g1 = repmat({'No Drug'},length(NDacc),1);
    g2 = repmat({'Bicuculine'},length(BICacc),1);
    g3 = repmat({'Muscimol'},length(MUSacc),1);
    g = [g1; g2; g3]; %Grouping
    xNetAcc = [NDacc BICacc MUSacc];
    mu_Net = mean(NDacc);
    mu_Bic = mean(BICacc);
    mu_Mus = mean(MUSacc);
    
end %End num TimeWin

%Figure out the smallest p values
[Pval_0Spon,indT_0Spon]=min(p_NDbic_0Spon);
[Wval_0Spon,indW_0Spon]=min(pW_NDbic_0Spon);
Tmins=[Pval_0Spon];
Wmins=[Wval_0Spon];
[ValSponTmin,SponTmin]=min(Tmins);
[ValSponWmin,SponWmin]=min(Wmins);
%Plot everything
ttestTtl=sprintf('Evoked Twin 2 Sample T-Test P-Values %s',odorName);
WCMttl=sprintf('Evoked Twin Wilcoxon Rank Sum Test P-Values %s',odorName);
ttestFile=sprintf('Ttest_Pvals_%s',odorName);
WCMfile=sprintf('WCMTest_Pvals_%s',odorName);
%Plot accuracy
f1=figure;
hold on
plot(p_NDbic_0Spon,'ko')
plot(p_NDmus_0Spon,'k*')
plot(0.01*ones(10,1),'k--')
xticklabels({'100ms','200ms','300ms','400ms','500ms','600ms','700ms','800ms','900ms','1000ms'})
legend('ND/Bic (0s)','ND/Mus (0s)','alpha = 0.01')
title(ttestTtl)
hold off
f2=figure;
hold on
plot(pW_NDbic_0Spon,'ko')
plot(pW_NDmus_0Spon,'k*')
plot(0.01*ones(10,1),'k--')
xticklabels({'100ms','200ms','300ms','400ms','500ms','600ms','700ms','800ms','900ms','1000ms'})
legend('ND/Bic (0s)','ND/Mus (0s)','alpha = 0.01')
title(WCMttl)
hold off

%% repeat but now remove first 300ms windows

fileName=sprintf('IndCellRAW3shift_OptWin_%s_%dsSpon.mat',odorName,sponTime);

p_NDbic_0Spon=[];
p_NDmus_0Spon=[];
pW_NDbic_0Spon=[];
pW_NDmus_0Spon=[];

for j=1:7
    PlotTtl=sprintf('All rat net accuracies - %s %dms',odorName,j*100);
    NDacc=[];
    BICacc=[];
    MUSacc=[];
    load(fileName,'optThrsh_ND','optThrsh_Bic','optThrsh_Mus')
    NDnumRats=size(optThrsh_ND,1);
    BICnumRats=size(optThrsh_Bic,1);
    MUSnumRats=size(optThrsh_Mus,1);
    for i=1:NDnumRats
        NDall=optThrsh_ND{i,4};
        NDacc=[NDacc NDall(j,:)];
    end
    for i=1:BICnumRats
        BICall=optThrsh_Bic{i,4};
        BICacc=[BICacc BICall(j,:)];
    end
    for i=1:MUSnumRats
        MUSall=optThrsh_Mus{i,4};
        MUSacc=[MUSacc MUSall(j,:)];
    end
    [~,p_NDbic]=ttest2(NDacc,BICacc,'VarType','unequal');
    [~,p_NDmus]=ttest2(NDacc,MUSacc,'VarType','unequal');
    pW_NDbic=ranksum(NDacc,BICacc);
    pW_NDmus=ranksum(NDacc,MUSacc);
    
    p_NDbic_0Spon=[p_NDbic_0Spon p_NDbic];
    p_NDmus_0Spon=[p_NDmus_0Spon p_NDmus];
    pW_NDbic_0Spon=[pW_NDbic_0Spon pW_NDbic];
    pW_NDmus_0Spon=[pW_NDmus_0Spon pW_NDmus];

    g1 = repmat({'No Drug'},length(NDacc),1);
    g2 = repmat({'Bicuculine'},length(BICacc),1);
    g3 = repmat({'Muscimol'},length(MUSacc),1);
    g = [g1; g2; g3]; %Grouping
    xNetAcc = [NDacc BICacc MUSacc];
    mu_Net = mean(NDacc);
    mu_Bic = mean(BICacc);
    mu_Mus = mean(MUSacc);
end %End num TimeWin

%Figure out the smallest p values
[Pval_0Spon,indT_0Spon]=min(p_NDbic_0Spon);
[Wval_0Spon,indW_0Spon]=min(pW_NDbic_0Spon);
Tmins=[Pval_0Spon];
Wmins=[Wval_0Spon];
[ValSponTmin,SponTmin]=min(Tmins);
[ValSponWmin,SponWmin]=min(Wmins);
%Plot everything
ttestTtl=sprintf('Evoked Twin 2 Sample T-Test P-Values %s',odorName);
WCMttl=sprintf('Evoked Twin Wilcoxon Rank Sum Test P-Values %s',odorName);
ttestFile=sprintf('Ttest_Pvals_%s',odorName);
WCMfile=sprintf('WCMTest_Pvals_%s',odorName);
%Plot accuracy
f1=figure;
hold on
plot(p_NDbic_0Spon,'ko')
plot(p_NDmus_0Spon,'k*')
plot(0.01*ones(7,1),'k--')
xticklabels({'400ms','500ms','600ms','700ms','800ms','900ms','1000ms'})
legend('ND/Bic (0s)','ND/Mus (0s)','alpha = 0.01')
title(ttestTtl)
hold off
f2=figure;
hold on
plot(pW_NDbic_0Spon,'ko')
plot(pW_NDmus_0Spon,'k*')
plot(0.01*ones(7,1),'k--')
xticklabels({'400ms','500ms','600ms','700ms','800ms','900ms','1000ms'})
legend('ND/Bic (0s)','ND/Mus (0s)','alpha = 0.01')
title(WCMttl)
hold off