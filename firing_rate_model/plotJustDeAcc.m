% plot decoding accuracies from dJustDeAcc.mat from model

load dJustDeAcc.mat

lneO=length(epsO);
lneR=length(epsR);

load('dFitMeans.mat','Am_O','Am_R','dt')

nuOint=sum(Am_O,2)*dt/(900); %hard-code 900 to get Hz
nuRint=sum(Am_R,2)*dt/(900);

Ormn=nuOint+repmat(epsO',3,1);
Rtmn=nuRint+repmat(epsR',3,1);

OrMat1=repmat(Ormn(1,:)',lneR,1);
RtMat1=reshape(repmat(Rtmn(1,:),lneO,1),lneO*lneR,1);
distVc1=abs(RtMat1-OrMat1);
DaVc1=reshape(DeAcc(:,:,1),lneO*lneR,1);

OrMat2=repmat(Ormn(2,:)',lneR,1);
RtMat2=reshape(repmat(Rtmn(2,:),lneO,1),lneO*lneR,1);
distVc2=abs(RtMat2-OrMat2);
DaVc2=reshape(DeAcc(:,:,2),lneO*lneR,1);

OrMat3=repmat(Ormn(3,:)',lneR,1);
RtMat3=reshape(repmat(Rtmn(3,:),lneO,1),lneO*lneR,1);
distVc3=abs(RtMat3-OrMat3);
DaVc3=reshape(DeAcc(:,:,3),lneO*lneR,1);

figure
hold on
plot(distVc1,DaVc1,'k.','MarkerSize',22)
plot(distVc2,DaVc2,'.','color',[.5490 0 1],'MarkerSize',22)
plot(distVc3,DaVc3,'.','color',[0 .7843 0],'MarkerSize',22)

%find where pop means started at
irMn=find(epsR==0);
ioMn=find(epsO==0);

idBase=sub2ind([lneO lneR],ioMn,irMn); %get index for baseline

plot(distVc1(idBase),DaVc1(idBase),'ro','LineWidth',2)
plot(distVc2(idBase),DaVc2(idBase),'ro','LineWidth',2)
plot(distVc3(idBase),DaVc3(idBase),'ro','LineWidth',2)

set(gca,'FontSize',18)
xlabel('Distance between Means')
ylabel('Decoding Accuracy')