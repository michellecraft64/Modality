% plot decoding accuracies from dJustDeAcc.mat from model

load dJustDeAcc.mat

lneO=length(epsO);
lneR=length(epsR);

load('dFitMeans.mat','Am_O','Am_R','dt','aspPv')

nuOint=sum(Am_O,2)*dt/(900); %hard-code 900 to get Hz
nuRint=sum(Am_R,2)*dt/(900);

Ormn=nuOint+repmat(epsO',3,1);
Rtmn=nuRint+repmat(epsR',3,1);

OrMat1=repmat(Ormn(1,:)',lneR,1);
RtMat1=reshape(repmat(Rtmn(1,:),lneO,1),lneO*lneR,1);
ScalDistVc1=abs(sqrt(RtMat1*aspPv(1,2))-sqrt(OrMat1*aspPv(1,1))); %mean/std = sqrt(mu*rho)
DaVc1=reshape(DeAcc(:,:,1),lneO*lneR,1);

OrMat2=repmat(Ormn(2,:)',lneR,1);
RtMat2=reshape(repmat(Rtmn(2,:),lneO,1),lneO*lneR,1);
ScalDistVc2=abs(sqrt(RtMat2*aspPv(2,2))-sqrt(OrMat2*aspPv(2,1))); 
DaVc2=reshape(DeAcc(:,:,2),lneO*lneR,1);

OrMat3=repmat(Ormn(3,:)',lneR,1);
RtMat3=reshape(repmat(Rtmn(3,:),lneO,1),lneO*lneR,1);
ScalDistVc3=abs(sqrt(RtMat3*aspPv(3,2))-sqrt(OrMat3*aspPv(3,1)));
DaVc3=reshape(DeAcc(:,:,3),lneO*lneR,1);


%% scaled dist mean
figure
hold on
plot(ScalDistVc1,DaVc1,'k.','MarkerSize',22)
plot(ScalDistVc2,DaVc2,'.','color',[.5490 0 1],'MarkerSize',22)
plot(ScalDistVc3,DaVc3,'.','color',[0 .7843 0],'MarkerSize',22)

%find where pop means started at
irMn=find(epsR==0);
ioMn=find(epsO==0);

idBase=sub2ind([lneO lneR],ioMn,irMn); %get index for baseline

plot(ScalDistVc1(idBase),DaVc1(idBase),'ro','LineWidth',2)
plot(ScalDistVc2(idBase),DaVc2(idBase),'ro','LineWidth',2)
plot(ScalDistVc3(idBase),DaVc3(idBase),'ro','LineWidth',2)

set(gca,'FontSize',18)
xlabel('Scaled Distance between Means')
ylabel('Decoding Accuracy')

% decod accur histograms
ccD=[128.01 128.01 128.01; 0 200 0; 140 0 255]./255;
edg=(.5-.025:.05:1.025)';
bw=0.1;

figure('Renderer', 'Painters');
hold on
histogram(DaVc1,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(1,:),'LineStyle','none')
histogram(DaVc3,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(2,:),'LineStyle','none')
histogram(DaVc2,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor','none','EdgeColor',ccD(3,:),'LineWidth',2)
plot(mean(DaVc1),.28,'.','color',ccD(1,:),'MarkerSize',22)
plot(mean(DaVc3),.26,'.','color',ccD(2,:),'MarkerSize',22)
plot(mean(DaVc2),.27,'.','color',ccD(3,:),'MarkerSize',22)
box off
set(gca,'FontSize',18)
axis([.5 1 0 .28])

%showing correlations in figure: |muR-muO| vs Decode Accur
drgPrep={'ND';'Bic';'Mus'}; 
pearsonsCorr=[corr(ScalDistVc1,DaVc1) ; corr(ScalDistVc3,DaVc3); corr(ScalDistVc2,DaVc2)  ];
RankCorr=[corr(ScalDistVc1,DaVc1,'type','Spearman') ; corr(ScalDistVc3,DaVc3,'type','Spearman')  ; corr(ScalDistVc2,DaVc2,'type','Spearman') ];
T=table(drgPrep,pearsonsCorr,RankCorr)


