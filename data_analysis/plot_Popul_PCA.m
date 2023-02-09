%script to process/show results of get_popDecod_PCA.m, saved in PopDecod_pcaLDA.mat

% also plots pop size 

load PopDecod_pcaLDA.mat

% calc percent of var explained

prcVr_expl_ND=zeros(length(vrPC_ND),1);
prcVr_expl_Bic=zeros(length(vrPC_Bic),1); 
prcVr_expl_Mus=zeros(length(vrPC_Mus),1); 
for j=1:size(prcVr_expl_ND,1)
    prcVr_expl_ND(j,:)=sum(vrPC_ND{j}(1:2))./sum(vrPC_ND{j});
end
for j=1:size(prcVr_expl_Bic,1)
    prcVr_expl_Bic(j,:)=sum(vrPC_Bic{j}(1:2))./sum(vrPC_Bic{j});
end
for j=1:size(prcVr_expl_Mus,1)
    prcVr_expl_Mus(j,:)=sum(vrPC_Mus{j}(1:2))./sum(vrPC_Mus{j});
end

ccD=[128.01 128.01 128.01; 0 200 0; 140 0 255]./255;
edg=(.5-.025:.05:1.025)';
bw=0.1;

figure('Renderer', 'Painters');
hold on
histogram(DeAc2_ND,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(1,:),'LineStyle','none')
histogram(DeAc2_Bic,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(2,:),'LineStyle','none')
histogram(DeAc2_Mus,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor','none','EdgeColor',ccD(3,:),'LineWidth',2)
plot(mean(DeAc2_ND),.55,'.','color',ccD(1,:),'MarkerSize',22)
plot(mean(DeAc2_Bic),.55,'.','color',ccD(2,:),'MarkerSize',22)
plot(mean(DeAc2_Mus),.55,'.','color',ccD(3,:),'MarkerSize',22)
box off
set(gca,'FontSize',18)
set(gca,'XLim',[.5 1])


%Significance
%Two-sample T-test
[~,pT_NDbic]=ttest2(DeAc2_ND,DeAc2_Bic,'VarType','unequal');
[~,pT_NDmus]=ttest2(DeAc2_ND,DeAc2_Mus,'VarType','unequal');
[~,pT_BICmus]=ttest2(DeAc2_Bic,DeAc2_Mus,'VarType','unequal');
%WCM rank sum
pW_NDbic=ranksum(DeAc2_ND,DeAc2_Bic);
pW_NDmus=ranksum(DeAc2_ND,DeAc2_Mus);
pW_BICmus=ranksum(DeAc2_Bic,DeAc2_Mus);
%One-way ANOVA
g_NDbic=[ones(length(DeAc2_ND),1);2*ones(length(DeAc2_Bic),1)];
g_NDmus=[ones(length(DeAc2_ND),1);3*ones(length(DeAc2_Mus),1)];
g_BICmus=[ones(length(DeAc2_Bic),1);3*ones(length(DeAc2_Mus),1)];
pA_NDbic = anova1([DeAc2_ND; DeAc2_Bic],g_NDbic,'off');
pA_NDmus = anova1([DeAc2_ND; DeAc2_Mus],g_NDmus,'off');
pA_BICmus = anova1([DeAc2_Bic; DeAc2_Mus],g_BICmus,'off');



%% prcnt var explained 
figure
subplot(1,3,1)
hold on
plot(prcVr_expl_ND(:,1),'ko-')
set(gca,'FontSize',18)
ylabel('Percent Var Explained')
subplot(1,3,2)
hold on
plot(prcVr_expl_Bic(:,1),'ko-')
set(gca,'FontSize',18)
ylabel('Percent Var Explained')
subplot(1,3,3)
hold on
plot(prcVr_expl_Mus(:,1),'ko-')
set(gca,'FontSize',18)
ylabel('Percent Var Explained')

figure
hold on
plot(sort(prcVr_expl_ND(:,1)),'.-','color',ccD(1,:),'MarkerSize',32)
plot(sort(prcVr_expl_Bic(:,1)),'.-','color',ccD(2,:),'MarkerSize',32)
plot(sort(prcVr_expl_Mus(:,1)),'.-','color',ccD(3,:),'MarkerSize',32)
set(gca,'FontSize',18)
xlabel('Recording Index (sorted)')
ylabel('Percent Variance Explained by PCA')
set(gca,'XLim',[1-eps 29+eps])

%% is prc var explain related to decoding accur?
figure
subplot(1,3,1)
hold on
plot(prcVr_expl_ND(:,1),DeAc2_ND,'ko')
set(gca,'FontSize',18)
xlabel('Percent Var Explained (2 or 3)'); ylabel('Decoding Accuracy')
m=corrcoef([prcVr_expl_ND(:,1) DeAc2_ND]); m(1,2)
subplot(1,3,2)
hold on
plot(prcVr_expl_Bic(:,1),DeAc2_Bic,'ko')
set(gca,'FontSize',18)
xlabel('Percent Var Explained (2 or 3)'); ylabel('Decoding Accuracy')
m=corrcoef([prcVr_expl_Bic(:,1) DeAc2_Bic]); m(1,2)
subplot(1,3,3)
hold on
plot(prcVr_expl_Mus(:,1),DeAc2_Mus,'ko')
set(gca,'FontSize',18)
xlabel('Percent Var Explained (2 or 3)'); ylabel('Decoding Accuracy')
m=corrcoef([prcVr_expl_Mus(:,1) DeAc2_Mus]); m(1,2)

%% plot size of record for BOTH SVM, PCA
load dSizeCells_perRecord.mat
sz_ND=[];
sz_Bic=[];
sz_Mus=[];
for j=1:size(nmCells_ND,1)
    sz_ND=[sz_ND; nmCells_ND{j}];
end
for j=1:size(nmCells_Bic,1)
    sz_Bic=[sz_Bic; nmCells_Bic{j}];
end
for j=1:size(nmCells_Mus,1)
    sz_Mus=[sz_Mus; nmCells_Mus{j}];
end
figure
hold on
plot(sort(sz_ND),'.-','color',ccD(1,:),'MarkerSize',32)
plot(sort(sz_Bic),'.-','color',ccD(2,:),'MarkerSize',32)
plot(sort(sz_Mus),'.-','color',ccD(3,:),'MarkerSize',32)
set(gca,'FontSize',18)
xlabel('Recording Index (sorted)')
ylabel('Population Size')
set(gca,'XLim',[1-eps 29+eps])    

