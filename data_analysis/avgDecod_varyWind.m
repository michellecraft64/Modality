% script to plot avg decoding accuracies as window size varies
clear
load IndCellRAW_OptWin_EB_0sSpon.mat

numWind=10;

mn_ND=zeros(numWind,1);  std_ND=zeros(numWind,1); 
mn_Bic=zeros(numWind,1);  std_Bic=zeros(numWind,1); 
mn_Mus=zeros(numWind,1);  std_Mus=zeros(numWind,1); 


for indWind=1:numWind %index from 1,..,10; using 9=900ms evoked b/c largest p-val

    tmp1=[]; tmp2=[]; tmp3=[];

    for j=1:8 %total # rats, for ND
        tmp1=[tmp1; optThrsh_ND{j,4}(indWind,:).'];
        if j<=3 % bic & mus
            tmp2=[tmp2; optThrsh_Bic{j,4}(indWind,:).'];
            tmp3=[tmp3; optThrsh_Mus{j,4}(indWind,:).'];
        end
        if j==4 % bic
            tmp2=[tmp2; optThrsh_Bic{j,4}(indWind,:).'];
        end
    end
    
    mn_ND(indWind)=mean(tmp1); std_ND(indWind)=std(tmp1);
    mn_Bic(indWind)=mean(tmp2); std_Bic(indWind)=std(tmp2);
    mn_Mus(indWind)=mean(tmp3); std_Mus(indWind)=std(tmp3);
end

ccD=[128.01 128.01 128.01; 0 200 0; 140 0 255]./255;
stdFac=1;
figure
hold on
plot(1:10,mn_ND,'color',ccD(1,:),'LineWidth',2)
plot(1:10,mn_ND+stdFac*std_ND,'-','color',ccD(1,:))
plot(1:10,mn_ND-stdFac*std_ND,'-','color',ccD(1,:))
plot(1:10,mn_Bic,'color',ccD(2,:),'LineWidth',2)
plot(1:10,mn_Bic+stdFac*std_Bic,'-','color',ccD(2,:))
plot(1:10,mn_Bic-stdFac*std_Bic,'-','color',ccD(2,:))
plot(1:10,mn_Mus,'color',ccD(3,:),'LineWidth',2)
plot(1:10,mn_Mus+stdFac*std_Mus,'-','color',ccD(3,:))
plot(1:10,mn_Mus-stdFac*std_Mus,'-','color',ccD(3,:))
set(gca,'YLim',[.45 1])

%% repeat but for first 300ms removed
clear
load IndCellRAW3shift_OptWin_EB_0sSpon.mat

numWind=7; %only 400ms, to 1s

mn_ND=zeros(numWind,1);  std_ND=zeros(numWind,1); 
mn_Bic=zeros(numWind,1);  std_Bic=zeros(numWind,1); 
mn_Mus=zeros(numWind,1);  std_Mus=zeros(numWind,1); 


for indWind=1:numWind %index from 1,..,10; using 9=900ms evoked b/c largest p-val

    tmp1=[]; tmp2=[]; tmp3=[];

    for j=1:8 %total # rats, for ND
        tmp1=[tmp1; optThrsh_ND{j,4}(indWind,:).'];
        if j<=3 % bic & mus
            tmp2=[tmp2; optThrsh_Bic{j,4}(indWind,:).'];
            tmp3=[tmp3; optThrsh_Mus{j,4}(indWind,:).'];
        end
        if j==4 % bic
            tmp2=[tmp2; optThrsh_Bic{j,4}(indWind,:).'];
        end
    end
    
    mn_ND(indWind)=mean(tmp1); std_ND(indWind)=std(tmp1);
    mn_Bic(indWind)=mean(tmp2); std_Bic(indWind)=std(tmp2);
    mn_Mus(indWind)=mean(tmp3); std_Mus(indWind)=std(tmp3);
end

ccD=[128.01 128.01 128.01; 0 200 0; 140 0 255]./255;
stdFac=1;
figure
hold on
plot(4:10,mn_ND,'color',ccD(1,:),'LineWidth',2)
plot(4:10,mn_ND+stdFac*std_ND,'-','color',ccD(1,:))
plot(4:10,mn_ND-stdFac*std_ND,'-','color',ccD(1,:))
plot(4:10,mn_Bic,'color',ccD(2,:),'LineWidth',2)
plot(4:10,mn_Bic+stdFac*std_Bic,'-','color',ccD(2,:))
plot(4:10,mn_Bic-stdFac*std_Bic,'-','color',ccD(2,:))
plot(4:10,mn_Mus,'color',ccD(3,:),'LineWidth',2)
plot(4:10,mn_Mus+stdFac*std_Mus,'-','color',ccD(3,:))
plot(4:10,mn_Mus-stdFac*std_Mus,'-','color',ccD(3,:))
set(gca,'YLim',[.45 1])