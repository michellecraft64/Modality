%Script to calculate all types of significance on trial-averaged firing
%rate

%% Load Data %%
clear
indCellFile=sprintf('IndCell_rawFR_EB.mat');
load(indCellFile)
Zmn_Orth_ND=NormStruct{1,1}; %nOB x 10 trials
Zmn_Ret_ND=NormStruct{1,2};
Zmn_Orth_Bic=NormStruct{2,1};
Zmn_Ret_Bic=NormStruct{2,2};
Zmn_Orth_Mus=NormStruct{3,1};
Zmn_Ret_Mus=NormStruct{3,2};


% get Variances
vr_Zmn_Orth_ND=var(Zmn_Orth_ND,0,2);
vr_Zmn_Ret_ND=var(Zmn_Ret_ND,0,2);
vr_Zmn_Orth_Bic=var(Zmn_Orth_Bic,0,2);
vr_Zmn_Ret_Bic=var(Zmn_Ret_Bic,0,2);
vr_Zmn_Orth_Mus=var(Zmn_Orth_Mus,0,2);
vr_Zmn_Ret_Mus=var(Zmn_Ret_Mus,0,2);


% Two-Sample T-test %%
% Between Drugs
[~,p_vrOr_NDbic]=ttest2(vr_Zmn_Orth_ND,vr_Zmn_Orth_Bic,'Vartype','unequal');
[~,p_vrOr_BICmus]=ttest2(vr_Zmn_Orth_Mus,vr_Zmn_Orth_Bic,'Vartype','unequal');
[~,p_vrOr_NDmus]=ttest2(vr_Zmn_Orth_ND,vr_Zmn_Orth_Mus,'Vartype','unequal');
[~,p_vrRet_NDbic]=ttest2(vr_Zmn_Ret_ND,vr_Zmn_Ret_Bic,'Vartype','unequal');
[~,p_vrRet_BICmus]=ttest2(vr_Zmn_Ret_Mus,vr_Zmn_Ret_Bic,'Vartype','unequal');
[~,p_vrRet_NDmus]=ttest2(vr_Zmn_Ret_ND,vr_Zmn_Ret_Mus,'Vartype','unequal');

% Wilcoxon rank-sum test %%
% Between drugs
[pW_vrOr_NDbic,~,zW_o_NDbic]=ranksum(vr_Zmn_Orth_ND,vr_Zmn_Orth_Bic);
[pW_vrOr_BICmus,~,zW_o_BICmus]=ranksum(vr_Zmn_Orth_Mus,vr_Zmn_Orth_Bic);
[pW_vrOr_NDmus,~,zW_o_NDmus]=ranksum(vr_Zmn_Orth_ND,vr_Zmn_Orth_Mus);
[pW_vrRet_NDbic,~,zW_r_NDbic]=ranksum(vr_Zmn_Ret_ND,vr_Zmn_Ret_Bic);
[pW_vrRet_BICmus,~,zW_r_BICmus]=ranksum(vr_Zmn_Ret_Mus,vr_Zmn_Ret_Bic);
[pW_vrRet_NDmus,~,zW_r_NDmus]=ranksum(vr_Zmn_Ret_ND,vr_Zmn_Ret_Mus);

%One-way ANOVA, Between drugs only
    g_oNDBic=[ones(length(vr_Zmn_Orth_ND),1);2*ones(length(vr_Zmn_Orth_Bic),1)];
[pA_vrOr_NDBic,tblA_o_NDbic] = anova1([vr_Zmn_Orth_ND; vr_Zmn_Orth_Bic],g_oNDBic,'off');
    g_oMusBic=[ones(length(vr_Zmn_Orth_Mus),1);2*ones(length(vr_Zmn_Orth_Bic),1)];
[pA_vrOr_MusBic,tblA_o_BICmus] = anova1([vr_Zmn_Orth_Mus; vr_Zmn_Orth_Bic],g_oMusBic,'off');
    g_oNDMus=[ones(length(vr_Zmn_Orth_ND),1);2*ones(length(vr_Zmn_Orth_Mus),1)];
[pA_vrOr_NDMus,tblA_o_NDmus] = anova1([vr_Zmn_Orth_ND; vr_Zmn_Orth_Mus],g_oNDMus,'off');

    g_rNDBic=[ones(length(vr_Zmn_Ret_ND),1);2*ones(length(vr_Zmn_Ret_Bic),1)];
[pA_vrRet_NDBic,tblA_r_NDbic] = anova1([vr_Zmn_Ret_ND; vr_Zmn_Ret_Bic],g_rNDBic,'off');
    g_rMusBic=[ones(length(vr_Zmn_Ret_Mus),1);2*ones(length(vr_Zmn_Ret_Bic),1)];
[pA_vrRet_MusBic,tblA_r_BICmus] = anova1([vr_Zmn_Ret_Mus; vr_Zmn_Ret_Bic],g_rMusBic,'off');
    g_rNDMus=[ones(length(vr_Zmn_Ret_ND),1);2*ones(length(vr_Zmn_Ret_Mus),1)];
[pA_vrRet_NDMus,tblA_r_NDmus] = anova1([vr_Zmn_Ret_ND; vr_Zmn_Ret_Mus],g_rNDMus,'off');


%% table for raw mean vars
drgs={'ND';'Bic';'Mus'};
OrthTrialVar=[mean(vr_Zmn_Orth_ND); mean(vr_Zmn_Orth_Bic); mean(vr_Zmn_Orth_Mus)];
RetrTrialVar=[mean(vr_Zmn_Ret_ND); mean(vr_Zmn_Ret_Bic); mean(vr_Zmn_Ret_Mus)];
T_avgTrialVar=table(drgs,OrthTrialVar,RetrTrialVar)

%% tables to show info in manuscript
%showing correlations in figure: |muR-muO| vs Decode Accur

drgPrepComp={'Mus<ND';'Mus<Bic';'Bic<ND'}; 

%for ortho
ttestO=[p_vrOr_NDmus; p_vrOr_BICmus; p_vrOr_NDbic];
wilcRankO=[pW_vrOr_NDmus; pW_vrOr_BICmus; pW_vrOr_NDbic];
owanovaO=[pA_vrOr_NDMus; pA_vrOr_MusBic ; pA_vrOr_NDBic ];

T_Var_orth=table(drgPrepComp,ttestO,wilcRankO,owanovaO)

% effect sizes
numND=length(vr_Zmn_Orth_ND); numBic=length(vr_Zmn_Orth_Bic); numMus=length(vr_Zmn_Orth_Mus);
EffSiz_ttst=[ abs(mean(vr_Zmn_Orth_ND)-mean(vr_Zmn_Orth_Mus))/sqrt( ((numND-1)*var(vr_Zmn_Orth_ND)+(numMus-1)*var(vr_Zmn_Orth_Mus))/(numND+numMus-2) ) ; ...
    abs(mean(vr_Zmn_Orth_Bic)-mean(vr_Zmn_Orth_Mus))/sqrt( ((numBic-1)*var(vr_Zmn_Orth_Bic)+(numMus-1)*var(vr_Zmn_Orth_Mus))/(numBic+numMus-2) ); ...
    abs(mean(vr_Zmn_Orth_ND)-mean(vr_Zmn_Orth_Bic))/sqrt( ((numND-1)*var(vr_Zmn_Orth_ND)+(numBic-1)*var(vr_Zmn_Orth_Bic))/(numND+numBic-2) )];
EffSiz_wrst=abs([zW_o_NDmus.zval/sqrt(numND+numMus); zW_o_BICmus.zval/sqrt(numBic+numMus); zW_o_NDbic.zval/sqrt(numND+numBic)]);
EffSiz_owanova=[tblA_o_NDmus{6}/tblA_o_NDmus{8}; tblA_o_BICmus{6}/tblA_o_BICmus{8}; tblA_o_NDbic{6}/tblA_o_NDbic{8}];
T_EffSiz_Ort=table(drgPrepComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)


%for retro
ttestR=[p_vrRet_NDmus; p_vrRet_BICmus; p_vrRet_NDbic];
wilcRankR=[pW_vrRet_NDmus; pW_vrRet_BICmus; pW_vrRet_NDbic];
owanovaR=[pA_vrRet_NDMus; pA_vrRet_MusBic ; pA_vrRet_NDBic ];

T_Var_retr=table(drgPrepComp,ttestR,wilcRankR,owanovaR)

% effect sizes
EffSiz_ttst=[ abs(mean(vr_Zmn_Ret_ND)-mean(vr_Zmn_Ret_Mus))/sqrt( ((numND-1)*var(vr_Zmn_Ret_ND)+(numMus-1)*var(vr_Zmn_Ret_Mus))/(numND+numMus-2) ) ; ...
    abs(mean(vr_Zmn_Ret_Bic)-mean(vr_Zmn_Ret_Mus))/sqrt( ((numBic-1)*var(vr_Zmn_Ret_Bic)+(numMus-1)*var(vr_Zmn_Ret_Mus))/(numBic+numMus-2) ); ...
    abs(mean(vr_Zmn_Ret_ND)-mean(vr_Zmn_Ret_Bic))/sqrt( ((numND-1)*var(vr_Zmn_Ret_ND)+(numBic-1)*var(vr_Zmn_Ret_Bic))/(numND+numBic-2) )];
EffSiz_wrst=abs([zW_r_NDmus.zval/sqrt(numND+numMus); zW_r_BICmus.zval/sqrt(numBic+numMus); zW_r_NDbic.zval/sqrt(numND+numBic)]);
EffSiz_owanova=[tblA_r_NDmus{6}/tblA_r_NDmus{8}; tblA_r_BICmus{6}/tblA_r_BICmus{8}; tblA_r_NDbic{6}/tblA_r_NDbic{8}];
T_EffSiz_Ret=table(drgPrepComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)