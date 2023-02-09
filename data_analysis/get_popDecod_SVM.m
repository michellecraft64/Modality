%script to process pop coded SVM, using Bayesian approach to find best
%BoxConstraints & kernelScale (trick to get nonlinear decoding)

% Alive rats: 1, 2, 6, 7, 8, 9, 10, 11
% Rats for ND/Bic: 1, 2, 8, 11
% Rats for ND/Mus: 6, 9, 10

save_flag=0; %if this is =1, WILL OVERWRITE EXISTING PopDecod_Bayes_SVM.mat

numTrials=20;
theClass=[-1*ones(numTrials/2,1); ones(numTrials/2,1)]; %assume equal trials

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
        DeAc_obj_ND=zeros(29,1);
        DeAc_crss_ND=zeros(29,1);
        cellNums=nmCells_ND;
    case 2
        indRts=ind_Bic;
        flNameEnd='_IndCell_EB_Bic';
        % %OUTPUTS to save
        DeAc_obj_Bic=zeros(12,1);
        DeAc_crss_Bic=zeros(12,1);
        cellNums=nmCells_Bic;
    case 3
        indRts=ind_Mus;
        flNameEnd='_IndCell_EB_Mus';
        % %OUTPUTS to save
        DeAc_obj_Mus=zeros(12,1);
        DeAc_crss_Mus=zeros(12,1);
        cellNums=nmCells_Mus;
end

cnt=1; %keep track of total # DecodesAccs
cntR=1; %subs for i index, raw count
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

        %setup partition for cross-valid
        c = cvpartition(numTrials,'KFold',10);

        % Bayesian optim to find best hyperparam & minimal cross-valid loss
        opts = struct('CVPartition',c,'AcquisitionFunctionName','expected-improvement-plus');

        Mdl = fitcsvm(X,theClass,'KernelFunction','rbf', ...
            'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts);

        CVMdl1 = crossval(Mdl); %10-fold cross-validation

        switch DrugPrep
            case 1
                %cross-validated misclassification rate, theory from objec
                DeAc_obj_ND(cnt,1) = 1 - Mdl.HyperparameterOptimizationResults.MinEstimatedObjective;
                DeAc_crss_ND(cnt,1) = 1 - kfoldLoss(CVMdl1); %mis-classification rate
            case 2
                DeAc_obj_Bic(cnt,1) = 1 - Mdl.HyperparameterOptimizationResults.MinEstimatedObjective;
                DeAc_crss_Bic(cnt,1) = 1- kfoldLoss(CVMdl1); %mis-classification rate
            case 3
                DeAc_obj_Mus(cnt,1) = 1 - Mdl.HyperparameterOptimizationResults.MinEstimatedObjective;
                DeAc_crss_Mus(cnt,1) = 1 - kfoldLoss(CVMdl1); %mis-classification rate
        end

        pause(0.5)
        close all

        subCnt=subCnt+cellNums{cntR,1}(k); %updates to count
        cnt=cnt+1;
    end
    cntR=cntR+1;
end

if(save_flag==1)
    switch DrugPrep
        case 1
            save('PopDecod_Bayes_SVM.mat','DeAc_obj_ND','DeAc_crss_ND') %assuming first time that .mat exists
        case 2
            save('PopDecod_Bayes_SVM.mat','-append','DeAc_obj_Bic','DeAc_crss_Bic') %assuming .mat exists (not first time)
        case 3
            save('PopDecod_Bayes_SVM.mat','-append','DeAc_obj_Mus','DeAc_crss_Mus') %assuming .mat exists (not first time)
    end
end

end

