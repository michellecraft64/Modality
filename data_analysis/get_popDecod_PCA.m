%script to process pop coding via PCA, fixed number of dims

% Alive rats: 1, 2, 6, 7, 8, 9, 10, 11
% Rats for ND/Bic: 1, 2, 8, 11
% Rats for ND/Mus: 6, 9, 10

save_flag=01;  %if this is =1, WILL OVERWRITE EXISTING PopDecod_pcaLDA.mat

numTrials=20;
hlTrls=10;

% Define ND/Bic vs ND/Mus indexing
ind_ND=[1 2 6 7 8 9 10 11];
ind_Bic=[1 2 8 11];
ind_Mus=[6 9 10];

load dSizeCells_perRecord.mat %get size of each recording (segmented by rat #)

for DrugPrep=1:3 %1=ND, 2=Bic, 3=Mus

switch DrugPrep
    case 1
        indRts=ind_ND;
        flNameEnd='_IndCell_EB_NoDrug';
        % %OUTPUTS to save
        DeAc2_ND=zeros(29,1);
        cellNums=nmCells_ND;
        vrPC_ND=cell(29,1);
    case 2
        indRts=ind_Bic;
        flNameEnd='_IndCell_EB_Bic';
        % %OUTPUTS to save
        DeAc2_Bic=zeros(12,1);
        cellNums=nmCells_Bic;
        vrPC_Bic=cell(12,1);
    case 3
        indRts=ind_Mus;
        flNameEnd='_IndCell_EB_Mus';
        % %OUTPUTS to save
        DeAc2_Mus=zeros(12,1);
        cellNums=nmCells_Mus;
        vrPC_Mus=cell(12,1);
end

cnt=1; %keep track of total # DecodesAccs
cntR=1; %subs for i index, raw count

lbls=[zeros(hlTrls,1); ones(hlTrls,1)]; %0=Ortho, 1=Retro

for i=indRts
    

fileName=['Rat',num2str(i),flNameEnd,'.mat'];

%Load file, set parameters
load(fileName) %Ortho/Retro raw spike counts (lenTime, numTrials, nID)

%hard-code -11 to reflect first 900ms
tmEnd=TimeVars.LastEvok-11;

scOR=squeeze(sum(sOR(TimeVars.FirstEvok:tmEnd,:,:)));
scRT=squeeze(sum(sRET(TimeVars.FirstEvok:tmEnd,:,:)));


    subCnt=0;
    for k=1:length(cellNums{cntR,1}) %loop over diff recordings within the same rat
        X=[scOR(:,subCnt+1:subCnt+cellNums{cntR,1}(k)) ; scRT(:,subCnt+1:subCnt+cellNums{cntR,1}(k))];

        [~,scrs_,lat_j]=pca(X);
        scrs_OR=scrs_(1:hlTrls,:); scrs_RT=scrs_(hlTrls+1:numTrials,:); 
        
        MdLin2=fitcdiscr(scrs_(:,1:2),lbls);
        categ=predict(MdLin2,scrs_(:,1:2));
        decdAcc2=sum(categ==lbls)/numTrials;

        % save lat_j, explained variances

        switch DrugPrep
            case 1
                DeAc2_ND(cnt,1) = decdAcc2;
                vrPC_ND{cnt,1}=lat_j;
            case 2
                DeAc2_Bic(cnt,1) = decdAcc2;
                vrPC_Bic{cnt,1}=lat_j;
            case 3
                DeAc2_Mus(cnt,1) = decdAcc2;
                vrPC_Mus{cnt,1}=lat_j;
        end

        subCnt=subCnt+cellNums{cntR,1}(k); %updates to count
        cnt=cnt+1;
    end
    cntR=cntR+1;
end

if(save_flag==1)
    switch DrugPrep
        case 1
            save('PopDecod_pcaLDA.mat','DeAc2_ND','vrPC_ND')  %assuming first time that .mat exists
        case 2
            save('PopDecod_pcaLDA.mat','-append','DeAc2_Bic','vrPC_Bic') %assuming .mat exists (not first time)
        case 3
            save('PopDecod_pcaLDA.mat','-append','DeAc2_Mus','vrPC_Mus') %assuming .mat exists (not first time)
    end
end

end