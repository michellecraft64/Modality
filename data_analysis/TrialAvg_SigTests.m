%Script to calculate all types of significance on trial-averaged net spike
%rate
clear

%% Load Data %%
indCellFile=sprintf('IndCell_normFR_EB.mat');
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
pW_vrOr_NDbic=ranksum(vr_Zmn_Orth_ND,vr_Zmn_Orth_Bic);
pW_vrOr_BICmus=ranksum(vr_Zmn_Orth_Mus,vr_Zmn_Orth_Bic);
pW_vrOr_NDmus=ranksum(vr_Zmn_Orth_ND,vr_Zmn_Orth_Mus);
pW_vrRet_NDbic=ranksum(vr_Zmn_Ret_ND,vr_Zmn_Ret_Bic);
pW_vrRet_BICmus=ranksum(vr_Zmn_Ret_Mus,vr_Zmn_Ret_Bic);
pW_vrRet_NDmus=ranksum(vr_Zmn_Ret_ND,vr_Zmn_Ret_Mus);

%One-way ANOVA, Between drugs only
    g_oNDBic=[ones(length(vr_Zmn_Orth_ND),1);2*ones(length(vr_Zmn_Orth_Bic),1)];
pA_vrOr_NDBic = anova1([vr_Zmn_Orth_ND; vr_Zmn_Orth_Bic],g_oNDBic,'off');
    g_oMusBic=[ones(length(vr_Zmn_Orth_Mus),1);2*ones(length(vr_Zmn_Orth_Bic),1)];
pA_vrOr_MusBic = anova1([vr_Zmn_Orth_Mus; vr_Zmn_Orth_Bic],g_oMusBic,'off');
    g_oNDMus=[ones(length(vr_Zmn_Orth_ND),1);2*ones(length(vr_Zmn_Orth_Mus),1)];
pA_vrOr_NDMus = anova1([vr_Zmn_Orth_ND; vr_Zmn_Orth_Mus],g_oNDMus,'off');

    g_rNDBic=[ones(length(vr_Zmn_Ret_ND),1);2*ones(length(vr_Zmn_Ret_Bic),1)];
pA_vrRet_NDBic = anova1([vr_Zmn_Ret_ND; vr_Zmn_Ret_Bic],g_rNDBic,'off');
    g_rMusBic=[ones(length(vr_Zmn_Ret_Mus),1);2*ones(length(vr_Zmn_Ret_Bic),1)];
pA_vrRet_MusBic = anova1([vr_Zmn_Ret_Mus; vr_Zmn_Ret_Bic],g_rMusBic,'off');
    g_rNDMus=[ones(length(vr_Zmn_Ret_ND),1);2*ones(length(vr_Zmn_Ret_Mus),1)];
pA_vrRet_NDMus = anova1([vr_Zmn_Ret_ND; vr_Zmn_Ret_Mus],g_rNDMus,'off');



%% tables to show info in manuscript
%showing correlations in figure: |muR-muO| vs Decode Accur

drgPrepComp={'Mus<ND';'Mus<Bic';'Bic<ND'}; 

%for ortho
ttestO=[p_vrOr_NDmus; p_vrOr_BICmus; p_vrOr_NDbic];
wilcRankO=[pW_vrOr_NDmus; pW_vrOr_BICmus; pW_vrOr_NDbic];
owanovaO=[pA_vrOr_NDMus; pA_vrOr_MusBic ; pA_vrOr_NDBic ];

T_orth=table(drgPrepComp,ttestO,wilcRankO,owanovaO)

%for retro
ttestR=[p_vrRet_NDmus; p_vrRet_BICmus; p_vrRet_NDbic];
wilcRankR=[pW_vrRet_NDmus; pW_vrRet_BICmus; pW_vrRet_NDbic];
owanovaR=[pA_vrRet_NDMus; pA_vrRet_MusBic ; pA_vrRet_NDBic ];

T_retr=table(drgPrepComp,ttestR,wilcRankR,owanovaR)