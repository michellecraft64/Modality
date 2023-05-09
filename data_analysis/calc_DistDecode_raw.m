%script to plot  and calculate correl between |muR - muO| vs Decode Accur

load DistMeans_raw.mat 
load DecodAcc_raw.mat

%no drug
figure
semilogx(mnSpread_ND,DecA_ND,'o','color',[1 1 1]*.52)
axis([ 1e-2-eps 16.72 0.5 1])
box off

%Bic
figure
semilogx(mnSpread_Bic,DecA_Bic,'.','color',[0 .78125 0],'MarkerSize',24)
axis([ 1e-2-eps 8.96 0.5 1])
box off

%Mus
figure
semilogx(mnSpread_Mus,DecA_Mus,'o','color',[.45902 0 1])
axis([ 1e-2-eps 9.1477 0.5 1])
box off

drgPrep={'ND';'Bic';'Mus'}; 

%showing correlations in figure: |muR-muO| vs Decode Accur
pearsonsCorr=[corr(mnSpread_ND,DecA_ND) ; corr(mnSpread_Bic,DecA_Bic) ; corr(mnSpread_Mus,DecA_Mus)];
RankCorr=[corr(mnSpread_ND,DecA_ND,'type','Spearman') ; corr(mnSpread_Bic,DecA_Bic,'type','Spearman')  ; corr(mnSpread_Mus,DecA_Mus,'type','Spearman') ];
T=table(drgPrep,pearsonsCorr,RankCorr)

%% do again with mean/std
max_allX=500; %cap x-vals to mean
%no drug
x_ND=mnSpread_ND./stSpread_ND;
sum(x_ND>max_allX)+sum(isnan(x_ND))
x_ND(isnan(x_ND))=mnSpread_ND(isnan(x_ND));
x_ND(x_ND>max_allX)=mnSpread_ND(x_ND>max_allX);
figure
semilogx(x_ND,DecA_ND,'o','color',[1 1 1]*.52)
axis([ 1e-2-eps max_allX 0.5 1])
box off

%Bic
x_Bic=mnSpread_Bic./stSpread_Bic;
sum(x_Bic>max_allX)+sum(isnan(x_Bic))
x_Bic(isnan(x_Bic))=mnSpread_Bic(isnan(x_Bic));
x_Bic(x_Bic>max_allX)=mnSpread_Bic(x_Bic>max_allX);
figure
semilogx(x_Bic,DecA_Bic,'.','color',[0 .78125 0],'MarkerSize',24)
axis([ 1e-2-eps max_allX 0.5 1])
box off

%Mus
x_Mus=mnSpread_Mus./stSpread_Mus;
sum(x_Mus>max_allX)+sum(isnan(x_Mus)) %how many are 'bad'
x_Mus(isnan(x_Mus))=mnSpread_Mus(isnan(x_Mus));
x_Mus(x_Mus>max_allX)=mnSpread_Mus(x_Mus>max_allX);
figure
semilogx(x_Mus,DecA_Mus,'o','color',[.45902 0 1])
axis([ 1e-2-eps max_allX 0.5 1])
box off

drgPrep={'ND';'Bic';'Mus'}; 

%showing correlations in figure: |muR-muO| vs Decode Accur
pearsonsCorr=[corr(x_ND,DecA_ND) ; corr(x_Bic,DecA_Bic) ; corr(x_Mus,DecA_Mus)];
RankCorr=[corr(x_ND,DecA_ND,'type','Spearman') ; corr(x_Bic,DecA_Bic,'type','Spearman')  ; corr(x_Mus,DecA_Mus,'type','Spearman') ];
T=table(drgPrep,pearsonsCorr,RankCorr)