%script to plot  and calculate correl between |muR - muO| vs Decode Accur

load DistMeans.mat 
load DecodAcc.mat

%no drug
figure
semilogx(mnSpread_ND,DecA_ND,'o','color',[1 1 1]*.52)
axis([ 0.0003 16.72 0.5 1])

%Bic
figure
semilogx(mnSpread_Bic,DecA_Bic,'.','color',[0 .78125 0],'MarkerSize',24)
axis([ 0.0003 8.96 0.5 1])

%Mus
figure
semilogx(mnSpread_Mus,DecA_Mus,'o','color',[.45902 0 1])
axis([ 0.0828 9.1477 0.5 1])

drgPrep={'ND';'Bic';'Mus'}; 

%showing correlations in figure: |muR-muO| vs Decode Accur
pearsonsCorr=[corr(mnSpread_ND,DecA_ND) ; corr(mnSpread_Bic,DecA_Bic) ; corr(mnSpread_Mus,DecA_Mus)];
RankCorr=[corr(mnSpread_ND,DecA_ND,'type','Spearman') ; corr(mnSpread_Bic,DecA_Bic,'type','Spearman')  ; corr(mnSpread_Mus,DecA_Mus,'type','Spearman') ];
T=table(drgPrep,pearsonsCorr,RankCorr)
