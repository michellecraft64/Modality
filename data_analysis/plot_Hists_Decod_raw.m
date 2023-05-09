% script to plot histogram on top of each other, with RAW (not subtr spon fr)

%% plot total decoding accuracies; FIG 1C!
load DecodAcc_raw.mat

ccD=[128.01 128.01 128.01; 0 200 0; 140 0 255]./255;
edg=(.5-.025:.05:1.025)';
bw=0.1;

figure('Renderer', 'Painters');
hold on
histogram(DecA_ND,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(1,:),'LineStyle','none')
histogram(DecA_Bic,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(2,:),'LineStyle','none')
histogram(DecA_Mus,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor','none','EdgeColor',ccD(3,:),'LineWidth',2)
plot(mean(DecA_ND),.23,'.','color',ccD(1,:),'MarkerSize',22)
plot(mean(DecA_Bic),.22,'.','color',ccD(2,:),'MarkerSize',22)
plot(mean(DecA_Mus),.245,'.','color',ccD(3,:),'MarkerSize',22)
box off
set(gca,'FontSize',18)
set(gca,'XLim',[.5 1])

%% p-vals for significance
[~,pT_NDbic]=ttest2(DecA_ND,DecA_Bic,'VarType','unequal');
[~,pT_NDmus]=ttest2(DecA_ND,DecA_Mus,'VarType','unequal');
%WCM rank sum
[pW_NDbic,~,zW_NDbic]=ranksum(DecA_ND,DecA_Bic);
[pW_NDmus,~,zW_NDmus]=ranksum(DecA_ND,DecA_Mus);
%One-way ANOVA
g_NDbic=[ones(length(DecA_ND),1);2*ones(length(DecA_Bic),1)];
g_NDmus=[ones(length(DecA_ND),1);2*ones(length(DecA_Mus),1)];
[pA_NDbic,tblA_NDbic] = anova1([DecA_ND; DecA_Bic],g_NDbic,'off');
[pA_NDmus,tblA_NDmus] = anova1([DecA_ND; DecA_Mus],g_NDmus,'off');

% table to show ND is stat signif better than Bic/Mus
drgComp={'Bic<ND';'Mus<ND'}; 
Ttest=[pT_NDbic; pT_NDmus]; 
WRankSum=[pW_NDbic; pW_NDmus];
owanova=[pA_NDbic; pA_NDmus];
T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

numND=length(DecA_ND); numBic=length(DecA_Bic); numMus=length(DecA_Mus);

EffSiz_ttst=[ abs(mean(DecA_ND)-mean(DecA_Bic))/sqrt( ((numND-1)*var(DecA_ND)+(numBic-1)*var(DecA_Bic))/(numND+numBic-2) ) ; ...
    abs(mean(DecA_ND)-mean(DecA_Mus))/sqrt( ((numND-1)*var(DecA_ND)+(numMus-1)*var(DecA_Mus))/(numND+numMus-2) )];
EffSiz_wrst=[zW_NDbic.zval/sqrt(numND+numBic); zW_NDmus.zval/sqrt(numND+numMus)];
EffSiz_owanova=[tblA_NDbic{6}/tblA_NDbic{8}; tblA_NDmus{6}/tblA_NDmus{8}];
T_EffSiz_DA=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)

%% plot prediction; condition on when mu_O<mu_R & mu_R<mu_O
load DecAcc_mnGreater_raw.mat 
load DecodAcc_raw.mat
ccD=[128.01 128.01 128.01; 0 200 0; 140 0 255]./255;
ccOg=[0 0 0; 0 .5 0; .275 0 .5];
edg=(.5-.025:.05:1.025)';
bw=0.1;
figure('Renderer', 'Painters');
hold on
histogram(De_rGo_Bic,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(2,:),'LineStyle','none')
histogram(De_rGo_ND,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(1,:),'LineStyle','none')
plot(mean(De_rGo_Bic),.21,'.','color',ccD(2,:),'MarkerSize',22)
plot(mean(De_rGo_ND),.21,'.','color',ccD(1,:),'MarkerSize',22)
box off
set(gca,'FontSize',18)
set(gca,'XLim',[.5 1])

figure('Renderer', 'Painters');
hold on
histogram(De_oGr_Mus,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(3,:),'LineStyle','none')
histogram(De_oGr_ND,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(1,:),'LineStyle','none')
plot(mean(De_oGr_Mus),.28,'.','color',ccD(3,:),'MarkerSize',22)
plot(mean(De_oGr_ND),.3,'.','color',ccD(1,:),'MarkerSize',22)
box off
set(gca,'FontSize',18)
set(gca,'XLim',[.5 1])
 
% these dont hold us as well (bottom row)
figure('Renderer', 'Painters');
hold on
histogram(De_rGo_Mus,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(3,:),'LineStyle','none')
histogram(De_rGo_ND,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(1,:),'LineStyle','none')
histogram(DecA_Mus,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor','none','EdgeColor',ccOg(3,:),'LineWidth',2)
plot(mean(De_rGo_Mus),.25,'.','color',ccD(3,:),'MarkerSize',22)
plot(mean(De_rGo_ND),.27,'.','color',ccD(1,:),'MarkerSize',22)
plot(mean(DecA_Mus),.25,'.','color',ccOg(3,:),'MarkerSize',22)
box off
set(gca,'FontSize',18)
set(gca,'XLim',[.5 1])

figure('Renderer', 'Painters');
hold on
histogram(De_oGr_Bic,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(2,:),'LineStyle','none')
histogram(De_oGr_ND,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(1,:),'LineStyle','none')
plot(mean(De_oGr_Bic),.29,'.','color',ccD(2,:),'MarkerSize',22)
plot(mean(De_oGr_ND),.28,'.','color',ccD(1,:),'MarkerSize',22)
box off
set(gca,'FontSize',18)
set(gca,'XLim',[.5 1])

%Significance
%Two-sample T-test
[~,pT_oGr_NDbic]=ttest2(De_oGr_ND,De_oGr_Bic,'VarType','unequal');
[~,pT_rGo_NDbic]=ttest2(De_rGo_ND,De_rGo_Bic,'VarType','unequal');
[~,pT_oGr_NDmus]=ttest2(De_oGr_ND,De_oGr_Mus,'VarType','unequal');
[~,pT_rGo_NDmus]=ttest2(De_rGo_ND,De_rGo_Mus,'VarType','unequal');
%WCM rank sum
pW_oGr_NDbic=ranksum(De_oGr_ND,De_oGr_Bic);
pW_rGo_NDbic=ranksum(De_rGo_ND,De_rGo_Bic);
pW_oGr_NDmus=ranksum(De_oGr_ND,De_oGr_Mus);
pW_rGo_NDmus=ranksum(De_rGo_ND,De_rGo_Mus);
%One-way ANOVA
g_oGr_NDbic=[ones(length(De_oGr_ND),1);2*ones(length(De_oGr_Bic),1)];
g_rGo_NDbic=[ones(length(De_rGo_ND),1);2*ones(length(De_rGo_Bic),1)];
g_oGr_NDmus=[ones(length(De_oGr_ND),1);3*ones(length(De_oGr_Mus),1)];
g_rGo_NDmus=[ones(length(De_rGo_ND),1);3*ones(length(De_rGo_Mus),1)];
pA_oGr_NDbic = anova1([De_oGr_ND; De_oGr_Bic],g_oGr_NDbic,'off');
pA_rGo_NDbic = anova1([De_rGo_ND; De_rGo_Bic],g_rGo_NDbic,'off');
pA_oGr_NDmus = anova1([De_oGr_ND; De_oGr_Mus],g_oGr_NDmus,'off');
pA_rGo_NDmus = anova1([De_rGo_ND; De_rGo_Mus],g_rGo_NDmus,'off');

% Significance of all Mus with condition muR<muO
[~,pT_oGr_bothMus]=ttest2(DecA_Mus,De_oGr_Mus,'VarType','unequal');
pW_oGr_bothMus=ranksum(DecA_Mus,De_oGr_Mus);
g_oGr_musAll=[ones(length(DecA_Mus),1);2*ones(length(De_oGr_Mus),1)];
pA_oGr_bothMus=anova1([DecA_Mus; De_oGr_Mus],g_oGr_musAll,'off');


