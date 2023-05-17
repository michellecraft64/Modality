% script for testing idea of Poisson Spike counts as trial var 
% (WONT match data decoding accuracy)
% has 1 excit, 1 inhib; implem in WCnet_EIcells.m 

%% 
clear

%Select params
tdE_arr=10; %in ms
trE_arr=2; 
tw=100; %in ms
w_Parr=[0.2 0.25 0.15]; %0.2 original, 0.25 MUS, 0.15 BIC
crr=1; % crr param; turn off plasticity with crr=0
plot_Flag=1; %Flag to plot comparisons of parameter differences
PoisDAcc_Flag=1; %Flag to sim Poisson trial to trial differences

ParamMat=combvec([tdE_arr;trE_arr],w_Parr);
numParams=size(ParamMat,2);

%Time variables
TimeVars=struct('dt',0.01,'lenSpon',500,'lenEvok',2000,'lenTrn',1000); %opt Twin = 900ms evoked (2s spont)

% % weights below & lamO/lamR 
%with original lamO/lamR
alphMax_or_MC=80; %Constant multiplier for input to MC
alphMax_or_I=95; %Constant multiplier for input to I
alphMax_ret_MC=50; %Constant multiplier for input to MC
alphMax_ret_I=15; %Constant multiplier for input to I

%Define time vectors
dt=TimeVars.dt;
tFullEnd=TimeVars.lenSpon+TimeVars.lenEvok+TimeVars.lenTrn;
tEnd=TimeVars.lenSpon+TimeVars.lenEvok; %without transient time
t=(0:dt:tEnd)';
Lt=length(t); %same length as Tmev
trn=(0:dt:tFullEnd-tEnd-dt)';
lenTrn=length(trn); %transient time
tmeEvok=(0:dt:TimeVars.lenEvok)';
tmeSpon=(-TimeVars.lenSpon:dt:0)';
Tmev=[tmeSpon(1:end-1); tmeEvok]/1000;
lenFullTme=tFullEnd/dt+1;
%Ortho input
lamO=struct('tsft',0,'lmOsp',.14,'lmOevk',0.1035,'tauO',95,'tauO2',100,'ssV',.15); % based on Craft et al 21
%tmShift=1; %in ms
TimeVars.lenEvokS=tmeEvok;%[tmeEvok; tmeEvok(end)+(dt:dt:tmShift)']; %take nuevok out further by tmShift
%nmShift=length((dt:dt:tmShift)); %# elem must shift by
nuevokO=lamO.lmOevk*(-(TimeVars.lenEvokS+lamO.tsft).*exp(-(TimeVars.lenEvokS+lamO.tsft)/lamO.tauO) + (TimeVars.lenEvokS+lamO.tsft).*exp(-(TimeVars.lenEvokS+lamO.tsft)/lamO.tauO2))+lamO.ssV;
%nuevokO=nuevokO(nmShift+1:end); %shift over by tmShift
lam_O=[lamO.lmOsp*ones(length(tmeSpon)-1,1);nuevokO];
%Retro input
lamR=struct('tsft',200,'lmRspon',.14,'lmRevk',0.001,'tau_R',5,'tauR2',500,'ssV',.25); % based on Craft et al 21
nuevok_R=lamR.lmRevk*(-(TimeVars.lenEvokS+lamR.tsft).*exp(-(TimeVars.lenEvokS+lamR.tsft)/lamR.tau_R) + (TimeVars.lenEvokS+lamR.tsft).*exp(-(TimeVars.lenEvokS+lamR.tsft)/lamR.tauR2))+lamR.ssV;
idSplt=round(202.2/dt+1); %see where plot(Tmev,lam_O,Tmev,lam_R,'r') intersects
nuevok_R(1:idSplt)=nuevokO(1:idSplt); %split, first part of Retro like Ortho
lam_R=[lamR.lmRspon*ones(length(tmeSpon)-1,1);nuevok_R];
%plot(Tmev,lam_O,Tmev,lam_R,'r')

InputCurrent_Or=struct('Mc',lam_O*alphMax_or_MC,'Pg',lam_O*alphMax_or_I);
InputCurrent_Ret=struct('Mc',lam_R*alphMax_ret_MC,'Pg',lam_R*alphMax_ret_I);

xMin=struct('xM',0,'xP',0); %minimal current for firing

%Initialize output
Am_O=zeros(numParams,Lt);
Ap_O=zeros(numParams,Lt);
syEo=zeros(numParams,Lt);
syIo=zeros(numParams,Lt);
w_MP_O=zeros(numParams,Lt);
Am_R=zeros(numParams,Lt);
Ap_R=zeros(numParams,Lt);
syEr=zeros(numParams,Lt);
syIr=zeros(numParams,Lt);
w_MP_R=zeros(numParams,Lt);
% save transients too
Amt_O=zeros(numParams,lenTrn);
Apt_O=zeros(numParams,lenTrn);
w_MPt_O=zeros(numParams,lenTrn);
Amt_R=zeros(numParams,lenTrn);
Apt_R=zeros(numParams,lenTrn);
w_MPt_R=zeros(numParams,lenTrn);

%Iterate through params
for i=1:numParams
    %Define specific time constants and coupling strengths
    TimeConstants=struct('tau_M',10,'tau_P',5.5,'tdE',ParamMat(1,i),'trE',ParamMat(2,i),'tdI',10,'trI',2,'tw',tw);
    % td[E/I]=decay time of syn, tr[E/I]=rise time of syn
    CpldWx=struct('w_I',1,'w_P',ParamMat(3,i),'wmp0',1); 
    

    %call for Ortho
    [Am_O(i,:),Ap_O(i,:),syEo(i,:),syIo(i,:),w_MP_O(i,:),Amt_O(i,:),Apt_O(i,:),w_MPt_O(i,:)]=WC_EIcells(CpldWx,TimeVars,TimeConstants,xMin,InputCurrent_Or,crr);
    %call for Ret
    [Am_R(i,:),Ap_R(i,:),syEr(i,:),syIr(i,:),w_MP_R(i,:),Amt_R(i,:),Apt_R(i,:),w_MPt_R(i,:)]=WC_EIcells(CpldWx,TimeVars,TimeConstants,xMin,InputCurrent_Ret,crr);


end

%% Poisson Spike Counts for Trial Var
idZer=round((0-Tmev(1))/dt*1000+1); idEvD=round((0.9-Tmev(1))/dt*1000+1);
if(PoisDAcc_Flag==1)
    tic
    Ntr=50000; %aspirational p for neg binomial, p=1 is Poisson, p inverse prop to variance
    AccMat=zeros(numParams,4);
    for k=1:numParams
        
        nuO=squeeze(Am_O(k,idZer:idEvD));
        % Poisson
        meanOc=sum(nuO)*dt;  %COUNTS, not fr
        
        nuR=squeeze(Am_R(k,idZer:idEvD));
        meanRc=sum(nuR)*dt; 
        sumSpks_O_Close=poissrnd(meanOc,Ntr,1)./(1000*0.9); %convert counts to rates
        sumSpks_R_Close=poissrnd(meanRc,Ntr,1)./(1000*0.9); %convert counts to rates

        figure('Renderer', 'Painters');
        hold on
        histogram(sumSpks_O_Close,'Normalization','probability','FaceColor','b','EdgeColor','b','LineStyle','none')
        histogram(sumSpks_R_Close,'Normalization','probability','FaceColor','r','EdgeColor','r','LineStyle','none')
        box off
        %Set decode acc params
        allRespn=[sumSpks_O_Close;sumSpks_R_Close];
        mnOr=mean(sumSpks_O_Close,1);
        mnRet=mean(sumSpks_R_Close,1);
        %Output (trial-avg net responses (netRespn); all LDA accuracies (OrAcc/RetAcc); averaged LDA accuracies (NeuronClass))
        % ---- speed up code by only checking out sub-set of threshold ---
        if(mnRet>mnOr) 
            lwThr=prctile(sumSpks_O_Close,35); %35th percentile of smaller
            hiThr=prctile(sumSpks_R_Close,65); %65th percentile of larger
        else
            lwThr=prctile(sumSpks_R_Close,35); %35th percentile of smaller
            hiThr=prctile(sumSpks_O_Close,65); %65th percentile of larger
        end
        indT=(allRespn>=lwThr)&(allRespn<=hiThr);
        thrsh=allRespn(indT); %consider subset for threshold values
        Nthrs=length(thrsh);
        RetAcc=zeros(Nthrs,1);
        OrAcc=zeros(Nthrs,1);
        ttlAcc=zeros(Nthrs,1);
        for i=1:Nthrs %consider subset for threshold values
            Respn=ones(Ntr*2,1);
            if mnOr>mnRet
                Respn(allRespn>thrsh(i))=0;
            elseif mnRet>mnOr
                Respn(allRespn<thrsh(i))=0;
            else
                Respn(allRespn<thrsh(i))=0;
            end
            OrAcc(i)=1-sum(Respn(1:Ntr))/Ntr;
            RetAcc(i)=sum(Respn(Ntr+1:end))/Ntr;
            ttlAcc(i)=(OrAcc(i)+RetAcc(i))/2;
        end % ---- end of 'new code' ----
        [~,OptInd]=max(ttlAcc);
        AccMat(k,:)=[thrsh(OptInd) OrAcc(OptInd) RetAcc(OptInd) ttlAcc(OptInd)];
        % plot optim threshold
        my=gca;
        my=my.YLim(end);
        Yvc=(0: (my)/(100-1) : my)';
        plot(thrsh(OptInd)*ones(100,1),Yvc,'k','LineWidth',2)
        
    end %over all params
    AccMat(:,4) %show all decoding accuracies (total)
    toc
    save dFR_Poisson AccMat 
end

%% table in manuscript for simulated trial-var (model) = xxx

or_nd_Dep=sum(Am_O(1,idZer:idEvD))*dt/(1000*0.9);
or_bic_Dep=sum(Am_O(3,idZer:idEvD))*dt/(1000*0.9);
or_mus_Dep=sum(Am_O(2,idZer:idEvD))*dt/(1000*0.9);
    rt_nd_Dep=sum(Am_R(1,idZer:idEvD))*dt/(1000*0.9);
    rt_bic_Dep=sum(Am_R(3,idZer:idEvD))*dt/(1000*0.9);
    rt_mus_Dep=sum(Am_R(2,idZer:idEvD))*dt/(1000*0.9);

drgComp={'Ortho';'Retro'}; 
Muscimol=[or_mus_Dep; rt_mus_Dep]./(1000*.9); 
NoDrug=[or_nd_Dep; rt_nd_Dep]./(1000*.9);
Bicuculline=[or_bic_Dep; rt_bic_Dep]./(1000*.9);
T_trlVar=table(drgComp,Muscimol,NoDrug,Bicuculline)



