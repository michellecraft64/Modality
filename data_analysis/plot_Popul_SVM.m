%script to process/show results of get_popDecod_SVM.m, saved in PopDecod_Bayes_SVM.mat

load PopDecod_Bayes_SVM.mat

ccD=[128.01 128.01 128.01; 0 200 0; 140 0 255]./255;
edg=(.5-.025:.05:1.025)';
bw=0.1;

figure('Renderer', 'Painters');
hold on
histogram(DeAc_crss_ND,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(1,:),'LineStyle','none')
histogram(DeAc_crss_Bic,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(2,:),'LineStyle','none')
histogram(DeAc_crss_Mus,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor','none','EdgeColor',ccD(3,:),'LineWidth',2)
plot(mean(DeAc_crss_ND),.45,'.','color',ccD(1,:),'MarkerSize',22)
plot(mean(DeAc_crss_Bic),.45,'.','color',ccD(2,:),'MarkerSize',22)
plot(mean(DeAc_crss_Mus),.45,'.','color',ccD(3,:),'MarkerSize',22)
box off
set(gca,'FontSize',18)
set(gca,'XLim',[.5 1])


%Significance
%Two-sample T-test
[~,pT_NDbic]=ttest2(DeAc_crss_ND,DeAc_crss_Bic,'VarType','unequal');
[~,pT_NDmus]=ttest2(DeAc_crss_ND,DeAc_crss_Mus,'VarType','unequal');
[~,pT_BICmus]=ttest2(DeAc_crss_Bic,DeAc_crss_Mus,'VarType','unequal');
%WCM rank sum
pW_NDbic=ranksum(DeAc_crss_ND,DeAc_crss_Bic);
pW_NDmus=ranksum(DeAc_crss_ND,DeAc_crss_Mus);
pW_BICmus=ranksum(DeAc_crss_Bic,DeAc_crss_Mus);
%One-way ANOVA
g_NDbic=[ones(length(DeAc_crss_ND),1);2*ones(length(DeAc_crss_Bic),1)];
g_NDmus=[ones(length(DeAc_crss_ND),1);3*ones(length(DeAc_crss_Mus),1)];
g_BICmus=[ones(length(DeAc_crss_Bic),1);3*ones(length(DeAc_crss_Mus),1)];
pA_NDbic = anova1([DeAc_crss_ND; DeAc_crss_Bic],g_NDbic,'off');
pA_NDmus = anova1([DeAc_crss_ND; DeAc_crss_Mus],g_NDmus,'off');
pA_BICmus = anova1([DeAc_crss_Bic; DeAc_crss_Mus],g_BICmus,'off');


