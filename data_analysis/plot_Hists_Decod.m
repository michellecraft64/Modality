% script to plot histogram on top of each other, CL 6/15/22 
% matching manuscript figures in Feb 2023

%% plot total decoding accuracies; FIG 1C!
load DecodAcc.mat

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


%% plot prediction; condition on when mu_O<mu_R & mu_R<mu_O
load DecAcc_mnGreater.mat
load DecodAcc.mat
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

