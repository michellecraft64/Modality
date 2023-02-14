%Script to optimize window of time for decoding.
% Uses brute-force get optimal threshold value for each rat/neuron 
% pair to have optimal decoding accuracy of ortho/retro
% Scores are calculated using spike counts from various times in evoked 
% state subtracted by the mean spontaneous spike counts (entire 3s).

% Output structure consisting of:
% 1) Optimal thresholds for each rat/neuron
% 2) Ortho accuracy
% 3) Retro accuracy
% 4) Net accuracy (mean(orthoAcc,retroAcc))
% Saves the structure with the optimal time window where "optimal" for time
% window is defined as that which results in the max significant difference
% between drug preparations
clear
% Select parameters
tstWin=0:1:9; %Array of tested time windows (from 100ms to 1000ms)
odor_to_keep = 2; %ODOR: EB = 2 (Food odor); hexa = 1 (non-Food odor)
save_flag = 0; %Flag to save (=1) or not save (=0) workspace
%For generating labels
or='Ortho';
ret='Retro';
odorName='EB';

szOptWin=size(tstWin,2);
for drgState=1:3
switch drgState
    case 1
        drug_to_keep = 0;
        ind=[1 2 6 7 8 9 10 11]; %actual good rats (excluding dead rats 3,4,5)
        numRats=length(ind);
        drugName='NoDrug';
        %Output
        optThrsh_ND=cell(numRats,4); %For every rat, each cell includes individual neuron [Opt Threshold value, Or Acc, Ret Acc, Total Acc]
        NumCellsIDnet_ND=zeros(numRats,1); %Number of cells that had identical net accuracy (and used acc that lead to smallest diff btwn or/ret)
    case 2
        drug_to_keep = 1;
        ind=[1 2 8 11]; %actual good rats (excluding dead rats 4,5)
        numRats=length(ind);
        drugName='Bic';
        %Output
        optThrsh_Bic=cell(numRats,4); %For every rat, each cell includes individual neuron [Opt Threshold value, Or Acc, Ret Acc, Total Acc]
        NumCellsIDnet_Bic=zeros(numRats,1); %Number of cells that had identical net accuracy (and used acc that lead to smallest diff btwn or/ret)
    case 3
        drug_to_keep = 2;
        ind=[6 9 10]; %actual good rats (excluding dead rat 3)
        numRats=length(ind);
        drugName='Mus';
        %Output
        optThrsh_Mus=cell(numRats,4); %For every rat, each cell includes individual neuron [Opt Threshold value, Or Acc, Ret Acc, Total Acc]
        NumCellsIDnet_Mus=zeros(numRats,1); %Number of cells that had identical net accuracy (and used acc that lead to smallest diff btwn or/ret)
end
for ratNum=ind
    fileName=sprintf('Rat%d_IndCell_%s_%s.mat',ratNum,odorName,drugName);
    %Load file, set parameters
    load(fileName) %Ortho/Retro raw spike counts (lenTime, numTrials, nID)
    numTrials=size(sOR,2)+size(sRET,2); %per or/ret cell
    numEvok=TimeVars.numEvok;
    StimShift=TimeVars.StimShift;
    FirstEvok=TimeVars.FirstEvok;
    LastSpon=TimeVars.LastSpon;
    FirstSpon=LastSpon-9; %(1s) %LastSpon-19; %(2s) %TimeVars.FirstSpon; %(3s) !!!!HARD CODED spon!!
    nOB=size(sOR,3);
    if nOB~=size(sRET,3)
        keyboard
    end
    %Output (trial-avg net responses (netRespn); all LDA accuracies (OrAcc/RetAcc); averaged LDA accuracies (NeuronClass))
    RetAcc=zeros(szOptWin,nOB,numTrials);
    OrAcc=zeros(szOptWin,nOB,numTrials);
    ttlAcc=zeros(szOptWin,nOB,numTrials);
    thrsh=zeros(szOptWin,nOB,numTrials);
    switch drgState
        case 1 %No drug
            AccMat_ND=zeros(szOptWin,nOB,4);
        case 2 %Bicuculine
            AccMat_Bic=zeros(szOptWin,nOB,4);
        case 3 %Muscimol
            AccMat_Mus=zeros(szOptWin,nOB,4);
    end
    for j=1:szOptWin
        LastEvok=FirstEvok+tstWin(j); %TimeVars.LastEvok; %HARD CODED-NEEDS CHANGING IF TWIN CHANGED
        lenEvok=length(FirstEvok:1:LastEvok);
        lenSpon=length(FirstSpon:1:LastSpon);
        for k=1:nOB
            %Calculate evoked spike scores
            %Spontaneous mean subtract retro
            tmpEvkRet = sRET(FirstEvok:LastEvok,:,k);
            tmpSpnRet = sRET(FirstSpon:LastSpon,:,k);
            %Spontaneous mean subtract ortho
            tmpEvkOrth = sOR(FirstEvok:LastEvok,:,k);
            tmpSpnOrth = sOR(FirstSpon:LastSpon,:,k);
            Zmn_RetTemp=sum(tmpEvkRet,1)/lenEvok-sum(tmpSpnRet,1)/lenSpon;
            Zmn_OrthTemp=sum(tmpEvkOrth,1)/lenEvok-sum(tmpSpnOrth,1)/lenSpon;
            % first 10 rows are ortho, next 10 (11-20) are retro
            OrRespn=Zmn_OrthTemp';
            RetRespn=Zmn_RetTemp';
            allRespn=[OrRespn;RetRespn];
            mnOr=mean(OrRespn,1);
            mnRet=mean(RetRespn,1);
            for i=1:numTrials
                Respn=ones(numTrials,1);
                if allRespn(i)==max(allRespn)
                    thrsh(j,k,i)=allRespn(i)-eps;
                else
                    thrsh(j,k,i)=allRespn(i)+eps;
                end
                if mnOr>mnRet
                    Respn(allRespn>thrsh(j,k,i))=0;
                    upperTxt='Ortho';
                    lowerTxt='Retro';
                elseif mnRet>mnOr
                    Respn(allRespn<thrsh(j,k,i))=0;
                    lowerTxt='Ortho';
                    upperTxt='Retro';
                else
                    Respn(allRespn<thrsh(j,k,i))=0;
                    lowerTxt='Ortho';
                    upperTxt='Retro';
                end
                OrAcc(j,k,i)=1-sum(Respn(1:10))/10;
                RetAcc(j,k,i)=sum(Respn(11:20))/10;
                ttlAcc(j,k,i)=(OrAcc(j,k,i)+RetAcc(j,k,i))/2;
            end
            OptInd=find(ttlAcc(j,k,:)==max(ttlAcc(j,k,:)));
            if length(OptInd)>1
                switch drgState
                    case 1 %No Drug
                        NumCellsIDnet_ND(ind==ratNum)=NumCellsIDnet_ND(ind==ratNum)+1;
                    case 2 %Bicuculine
                        NumCellsIDnet_Bic(ind==ratNum)=NumCellsIDnet_Bic(ind==ratNum)+1;    
                    case 3 %Muscimol
                        NumCellsIDnet_Mus(ind==ratNum)=NumCellsIDnet_Mus(ind==ratNum)+1;
                end
                chkDiff=diff([OrAcc(j,k,OptInd);RetAcc(j,k,OptInd)]);
                [~,minInd]=min(chkDiff);
                OptInd=OptInd(minInd); 
            end
            switch drgState
                case 1 %No drugs
                    AccMat_ND(j,k,:)=[thrsh(j,k,OptInd) OrAcc(j,k,OptInd) RetAcc(j,k,OptInd) ttlAcc(j,k,OptInd)];
                case 2 %Bicuculine
                    AccMat_Bic(j,k,:)=[thrsh(j,k,OptInd) OrAcc(j,k,OptInd) RetAcc(j,k,OptInd) ttlAcc(j,k,OptInd)];
                case 3 %Muscimol
                    AccMat_Mus(j,k,:)=[thrsh(j,k,OptInd) OrAcc(j,k,OptInd) RetAcc(j,k,OptInd) ttlAcc(j,k,OptInd)];
            end
        end
    end
    switch drgState
        case 1 %No Drug
            optThrsh_ND{ratNum,1}=AccMat_ND(:,:,1);
            optThrsh_ND{ratNum,2}=AccMat_ND(:,:,2);
            optThrsh_ND{ratNum,3}=AccMat_ND(:,:,3);
            optThrsh_ND{ratNum,4}=AccMat_ND(:,:,4);
            optThrsh_ND=optThrsh_ND(~any(cellfun('isempty', optThrsh_ND), 2), :);
        case 2 %Bicuculine
            optThrsh_Bic{ratNum,1}=AccMat_Bic(:,:,1);
            optThrsh_Bic{ratNum,2}=AccMat_Bic(:,:,2);
            optThrsh_Bic{ratNum,3}=AccMat_Bic(:,:,3);
            optThrsh_Bic{ratNum,4}=AccMat_Bic(:,:,4);
            optThrsh_Bic=optThrsh_Bic(~any(cellfun('isempty', optThrsh_Bic), 2), :);
        case 3 %Muscimol
            optThrsh_Mus{ratNum,1}=AccMat_Mus(:,:,1);
            optThrsh_Mus{ratNum,2}=AccMat_Mus(:,:,2);
            optThrsh_Mus{ratNum,3}=AccMat_Mus(:,:,3);
            optThrsh_Mus{ratNum,4}=AccMat_Mus(:,:,4);
            optThrsh_Mus=optThrsh_Mus(~any(cellfun('isempty', optThrsh_Mus), 2), :);
    end
end %end loop all rats
end %end switch drug state
%Two-sample F-test for equal variances
[~,pVar_NDbic] = vartest2(AccMat_ND(:,:,4)',AccMat_Bic(:,:,4)');
[~,pVar_NDmus] = vartest2(AccMat_ND(:,:,4)',AccMat_Mus(:,:,4)');
%Two-sample T-test
[~,p_NDbic]=ttest2(AccMat_ND(:,:,4)',AccMat_Bic(:,:,4)','VarType','unequal');
[~,p_NDmus]=ttest2(AccMat_ND(:,:,4)',AccMat_Mus(:,:,4)','VarType','unequal');
%Wilcoxon-Mann-Whitney test (rank sum)
S(1).ND=AccMat_ND(1,:,4);S(2).ND=AccMat_ND(2,:,4);S(3).ND=AccMat_ND(3,:,4);
S(4).ND=AccMat_ND(4,:,4);S(5).ND=AccMat_ND(5,:,4);S(6).ND=AccMat_ND(6,:,4);
S(7).ND=AccMat_ND(7,:,4);S(8).ND=AccMat_ND(8,:,4);S(9).ND=AccMat_ND(9,:,4);
S(10).ND=AccMat_ND(10,:,4);
S(1).Bic=AccMat_Bic(1,:,4);S(2).Bic=AccMat_Bic(2,:,4);S(3).Bic=AccMat_Bic(3,:,4);
S(4).Bic=AccMat_Bic(4,:,4);S(5).Bic=AccMat_Bic(5,:,4);S(6).Bic=AccMat_Bic(6,:,4);
S(7).Bic=AccMat_Bic(7,:,4);S(8).Bic=AccMat_Bic(8,:,4);S(9).Bic=AccMat_Bic(9,:,4);
S(10).Bic=AccMat_Bic(10,:,4);
S(1).Mus=AccMat_Mus(1,:,4);S(2).Mus=AccMat_Mus(2,:,4);S(3).Mus=AccMat_Mus(3,:,4);
S(4).Mus=AccMat_Mus(4,:,4);S(5).Mus=AccMat_Mus(5,:,4);S(6).Mus=AccMat_Mus(6,:,4);
S(7).Mus=AccMat_Mus(7,:,4);S(8).Mus=AccMat_Mus(8,:,4);S(9).Mus=AccMat_Mus(9,:,4);
S(10).Mus=AccMat_Mus(10,:,4);
pW_NDbic = arrayfun(@(a) ranksum(a.ND,a.Bic),S);
pW_NDmus = arrayfun(@(a) ranksum(a.ND,a.Mus),S);
%One-way ANOVA
g_NDbic=[ones(size(AccMat_ND(1,:,4),2),1);2*ones(size(AccMat_Bic(1,:,4),2),1)];
S(1).g_NDbic=g_NDbic;S(2).g_NDbic=g_NDbic;S(3).g_NDbic=g_NDbic;S(4).g_NDbic=g_NDbic;
S(5).g_NDbic=g_NDbic;S(6).g_NDbic=g_NDbic;S(7).g_NDbic=g_NDbic;S(8).g_NDbic=g_NDbic;
S(9).g_NDbic=g_NDbic;S(10).g_NDbic=g_NDbic;
g_NDmus=[ones(size(AccMat_ND(1,:,4),2),1);3*ones(size(AccMat_Mus(1,:,4),2),1)];
S(1).g_NDmus=g_NDmus;S(2).g_NDmus=g_NDmus;S(3).g_NDmus=g_NDmus;S(4).g_NDmus=g_NDmus;
S(5).g_NDmus=g_NDmus;S(6).g_NDmus=g_NDmus;S(7).g_NDmus=g_NDmus;S(8).g_NDmus=g_NDmus;
S(9).g_NDmus=g_NDmus;S(10).g_NDmus=g_NDmus;
pANOVA_NDbic = arrayfun(@(a) anova1([a.ND a.Bic],a.g_NDbic,'off'),S);
pANOVA_NDmus = arrayfun(@(a) anova1([a.ND a.Mus],a.g_NDmus,'off'),S);
[min_NDbic,optWinInd_NDbic]=min(p_NDbic);
[min_NDmus,optWinInd_NDmus]=min(p_NDmus);
if min_NDbic<min_NDmus && p_NDmus(optWinInd_NDbic)<0.01
    optWin=optWinInd_NDbic*100; %Optimal time window in ms - needs changing if spike count Twin changes
elseif min_NDbic>min_NDmus && p_NDbic(optWinInd_NDmus)<0.01
    optWin=optWinInd_NDmus*100; %Optimal time window in ms - needs changing if spike count Twin changes
else
    sigInd_NDbic=find(p_NDbic<0.01);
    sigInd_NDmus=find(p_NDmus<0.01);
    optWin=sigInd_NDbic(ismember(sigInd_NDbic,sigInd_NDmus)); %Optimal time window in ms - needs changing if spike count Twin changes
    if size(optWin,2)>1
        optWin = optWin(1);
    end
end
%Include save function
if save_flag==1
    WkspcNm=sprintf('IndCellThrsh_OptWin_%s_%dsSpon.mat',odorName,lenSpon/10);
    save(WkspcNm,'optThrsh_ND','optThrsh_Bic','optThrsh_Mus','optWin',...
        'p_NDbic','p_NDmus','pVar_NDbic','pVar_NDmus','pW_NDbic','pW_NDmus',...
        'pANOVA_NDbic','pANOVA_NDmus')
end