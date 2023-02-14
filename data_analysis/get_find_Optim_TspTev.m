%Script to plot and find Optimal Twin (t*)
clear
% Select parameters
odor_to_keep = 2; 
plotInd_flag = 0; %Flag to show each time window plots (=1) or not show (=0) plots
%For generating labels
or='Ortho';
ret='Retro';
odorName='EB';

for sponTime=1:3
switch sponTime
    case 1
        fileName=sprintf('IndCellThrsh_OptWin_%s_%dsSpon.mat',odorName,sponTime);
        p_NDbic_1sSpon=[];
        p_NDmus_1sSpon=[];
        pW_NDbic_1sSpon=[];
        pW_NDmus_1sSpon=[];
    case 2
        fileName=sprintf('IndCellThrsh_OptWin_%s_%dsSpon.mat',odorName,sponTime);
        p_NDbic_2sSpon=[];
        p_NDmus_2sSpon=[];
        pW_NDbic_2sSpon=[];
        pW_NDmus_2sSpon=[];
    case 3
        fileName=sprintf('IndCellThrsh_OptWin_%s_%dsSpon.mat',odorName,sponTime);
        p_NDbic_3sSpon=[];
        p_NDmus_3sSpon=[];
        pW_NDbic_3sSpon=[];
        pW_NDmus_3sSpon=[];
end
for j=1:10
    PlotTtl=sprintf('All rat net accuracies - %s %dms %ds Spon',odorName,j*100,sponTime);
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
    switch sponTime
        case 1
            p_NDbic_1sSpon=[p_NDbic_1sSpon p_NDbic];
            p_NDmus_1sSpon=[p_NDmus_1sSpon p_NDmus];
            pW_NDbic_1sSpon=[pW_NDbic_1sSpon pW_NDbic];
            pW_NDmus_1sSpon=[pW_NDmus_1sSpon pW_NDmus];
        case 2
            p_NDbic_2sSpon=[p_NDbic_2sSpon p_NDbic];
            p_NDmus_2sSpon=[p_NDmus_2sSpon p_NDmus];
            pW_NDbic_2sSpon=[pW_NDbic_2sSpon pW_NDbic];
            pW_NDmus_2sSpon=[pW_NDmus_2sSpon pW_NDmus];
        case 3
            p_NDbic_3sSpon=[p_NDbic_3sSpon p_NDbic];
            p_NDmus_3sSpon=[p_NDmus_3sSpon p_NDmus];
            pW_NDbic_3sSpon=[pW_NDbic_3sSpon pW_NDbic];
            pW_NDmus_3sSpon=[pW_NDmus_3sSpon pW_NDmus];
    end
    g1 = repmat({'No Drug'},length(NDacc),1);
    g2 = repmat({'Bicuculine'},length(BICacc),1);
    g3 = repmat({'Muscimol'},length(MUSacc),1);
    g = [g1; g2; g3]; %Grouping
    xNetAcc = [NDacc BICacc MUSacc];
    mu_Net = mean(NDacc);
    mu_Bic = mean(BICacc);
    mu_Mus = mean(MUSacc);
    if plotInd_flag==1
        f0=figure;
        hold on
        boxplot(xNetAcc,g) %FIX THIS
        plot(1,mu_Net,'dg')
        text(1,mu_Net,sprintf('mu=%.2f',mu_Net))
        plot(2,mu_Bic,'dg')
        text(2,mu_Bic,sprintf('mu=%.2f',mu_Bic))
        plot(3,mu_Mus,'dg')
        text(3,mu_Mus,sprintf('mu=%.2f',mu_Mus))
        ylabel('Accuracies')
        title(PlotTtl)
        hold off
    end
end %End num TimeWin
end %End switch sponTime
%Figure out the smallest p values
[Pval_1sSpon,indT_1sSpon]=min(p_NDbic_1sSpon);
[Pval_2sSpon,indT_2sSpon]=min(p_NDbic_2sSpon);
[Pval_3sSpon,indT_3sSpon]=min(p_NDbic_3sSpon);
[Wval_1sSpon,indW_1sSpon]=min(pW_NDbic_1sSpon);
[Wval_2sSpon,indW_2sSpon]=min(pW_NDbic_2sSpon);
[Wval_3sSpon,indW_3sSpon]=min(pW_NDbic_3sSpon);
Tmins=[Pval_1sSpon Pval_2sSpon Pval_3sSpon];
Wmins=[Wval_1sSpon Wval_2sSpon Wval_3sSpon];
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
plot(p_NDbic_1sSpon,'go')
plot(p_NDmus_1sSpon,'g*')
plot(p_NDbic_2sSpon,'co')
plot(p_NDmus_2sSpon,'c*')
plot(p_NDbic_3sSpon,'mo')
plot(p_NDmus_2sSpon,'m*')
plot(0.01*ones(10,1),'k--')
xticklabels({'100ms','200ms','300ms','400ms','500ms','600ms','700ms','800ms','900ms','1000ms'})
legend('ND/Bic (1s)','ND/Mus (1s)','ND/Bic (2s)','ND/Mus (2s)','ND/Bic (3s)','ND/Mus (3s)','alpha = 0.01')
title(ttestTtl)
hold off
f2=figure;
hold on
plot(pW_NDbic_1sSpon,'go')
plot(pW_NDmus_1sSpon,'g*')
plot(pW_NDbic_2sSpon,'co')
plot(pW_NDmus_2sSpon,'c*')
plot(pW_NDbic_3sSpon,'mo')
plot(pW_NDmus_2sSpon,'m*')
plot(0.01*ones(10,1),'k--')
xticklabels({'100ms','200ms','300ms','400ms','500ms','600ms','700ms','800ms','900ms','1000ms'})
legend('ND/Bic (1s)','ND/Mus (1s)','ND/Bic (2s)','ND/Mus (2s)','ND/Bic (3s)','ND/Mus (3s)','alpha = 0.01')
title(WCMttl)
hold off
