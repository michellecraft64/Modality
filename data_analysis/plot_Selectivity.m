% plots fraction of orth / retro/ neither, all based on trial-avg spike  counts, for anesthetized
% for Fig 1Aiii

load dSizeCells_perRecord.mat %get size of each recording (segmented by rat #)

numTrials=20;
hlTrls=10;

numRecs=[29; 12; 12]; %number of recordings
nmCells=[913; 413; 419];

% Define ND/Bic vs ND/Mus indexing
ind_ND=[1 2 6 7 8 9 10 11];
ind_Bic=[1 2 8 11];
ind_Mus=[6 9 10];

% %OUTPUTS to save
scTrl=cell(3,1);

for dp_j=1:3 %1=ND, 2=Bic, 3=Mus
    scTrl{dp_j}=zeros(nmCells(dp_j),2); %1st col=orth, 2nd col=retr
    
    temp=[];

switch dp_j
    case 1
        indRts=ind_ND;
        flNameEnd='_IndCell_EB_NoDrug';
        cellNums=nmCells_ND;
    case 2
        indRts=ind_Bic;
        flNameEnd='_IndCell_EB_Bic';
        cellNums=nmCells_Bic;
    case 3
        indRts=ind_Mus;
        flNameEnd='_IndCell_EB_Mus';
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

        temp=[temp; mean(X(1:hlTrls,:))' mean(X(hlTrls+1:end,:))'];
        
        subCnt=subCnt+cellNums{cntR,1}(k); %updates to count
        cnt=cnt+1;
    end
    cntR=cntR+1;
end

scTrl{dp_j}=temp;

end

%% plots for paper
num_Prefe=zeros(3,3); %row=drug (1=ND,2=Bic,3=Mus), col=Prefe (1=Ortho, 2=Neither, 3=Retro, 
for j=1:3
    scn_srted=sort(scTrl{j},2); %sorted; small in (:,1), largest in (:,2)
    id_Neith= ( scn_srted(:,2)-scn_srted(:,1) <= (0.01*scn_srted(:,2)) ); %which MCs within 1%
    num_Prefe(j,2)=sum(id_Neith);
    %num_Prefe(j,3)/length(id_Neith)
    
    id_Orth=(scTrl{j}(:,1)>scTrl{j}(:,2)).*(id_Neith==0);  %get indices of all orth pref cells 
    num_Prefe(j,1)=sum(id_Orth);
    id_Retr=(scTrl{j}(:,2)>scTrl{j}(:,1)).*(id_Neith==0);  %get indices of all retr pref cells 
    num_Prefe(j,3)=sum(id_Retr);
end

cumPerct = cumsum(num_Prefe,2)./repmat(nmCells,1,3);
xx=(0:.1:1)';
lnX=length(xx);
figure
hold on
for j=1:3
    plot(xx,j*ones(lnX,1),'k')
    plot(xx,(j+.4)*ones(lnX,1),'k')
    plot(cumPerct(j,1)*ones(1,5), j:.1:j+.4,'r')
    plot(cumPerct(j,2)*ones(1,5), j:.1:j+.4,'r')
end
axis([-.2 1.2 1-.2 3.6])

