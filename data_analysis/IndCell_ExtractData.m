%script to create mat files Rat[j]_IndCell_EB_[NoDrug/Bic/Mus].mat
clear
basedir = pwd;
addpath([basedir '/code_data_analysis/']); 

SaveFlag=0; %Flag to save (=1) or not save (=0) .mat data files

odor_to_keep = 2;
odorName='EB';

prestr = '../experimental_data/';
% -- 11 rats total, only use 8 --
ratst_cell=cell(11,1);
ratst_cell{1}='011416';
ratst_cell{2}='012516'; 
ratst_cell{3}='020416';  %not as responsive, died, DONT USE
ratst_cell{4}='020916';  %not as responsive, died, DONT USE
ratst_cell{5}='021116';  %not as responsive, died, DONT USE
ratst_cell{6}='021816';
ratst_cell{7}='040115';
ratst_cell{8}='040515';
ratst_cell{9}='040915';
ratst_cell{10}='041615';
ratst_cell{11}='042115';

% Parameters you may want to change
Twin = 0.1; %this has to divide 2 evenly
SponL = 3;
Toverlap = 0; %MUST BE = 0 FOR THIS SCRIPT otherwise 0=disjoint, 1=hlfOverlapping

% Compute time-dependent
%     setParam:    Info for calc
%       setParam.Twin:            Time window used for PSTH
%       setParam.SponSampLength: Length of time before stim to use for
%               sampling spontaneous activity
%
setParam = [];
setParam.Twin       = Twin;
setParam.SponSampLength =SponL;
setParam.Toverlap=Toverlap;

for drgState=1:3
%drug_to_keep = [0]; %!! with drug (~=0), throw out 2 rats to be conservative!!
%  DRUG
%    no drug   = 0
%    Bicu      = 1  %blocks GABA_a inhib. synapses
%    Musc      = 2  %increases effects GABA_a
switch drgState
    case 1
        drug_to_keep = 0;
        ind=[1 2 6 7 8 9 10 11]; %all 8 good rats 
    case 2
        drug_to_keep = 1;
        ind=[1 2 8 11]; % 4 good rats
    case 3
        drug_to_keep = 2;
        ind=[6 9 10]; % 3 good rats
end
for ind_Rat=ind
    %Getting all files in directory
    datfilestr = dir([prestr ratst_cell{ind_Rat} '/' '*.mat']); %Getting all files in directory
    
    % Get metadata
    metaData = getMetaData_OdorDrug_SngleRat_fn(datfilestr,[prestr ratst_cell{ind_Rat}]);
    
    datfilelist = {};
    for j1=1:length(datfilestr)
        datfilelist = [datfilelist [prestr ratst_cell{ind_Rat} '/' datfilestr(j1).name]];
    end
    
    % Now restrict file list to ONLY keep desired files
    HasOdor = sum(metaData(:,2)==odor_to_keep,2);
    HasDrug = sum(metaData(:,3)==drug_to_keep,2);
    
    wantedFiles = find(HasOdor & HasDrug);
    
    datfilelist = datfilelist(wantedFiles);

    [SpikesOR,SpikesRET,TimeVars]=getFR_indcells(datfilelist,setParam);
    s=[];
    s.sOR=SpikesOR;
    s.sRET=SpikesRET;
    s.TimeVars=TimeVars;
    s.setParam=setParam;
    switch drgState
        case 1
            fileName=sprintf('Rat%d_IndCell_%s_NoDrug.mat',ind_Rat,odorName);
        case 2
            fileName=sprintf('Rat%d_IndCell_%s_Bic.mat',ind_Rat,odorName);
        case 3
            fileName=sprintf('Rat%d_IndCell_%s_Mus.mat',ind_Rat,odorName);
    end
    if(SaveFlag==1)
        save(fileName, '-struct', 's')
    end
end %for-loop over all rats

end %for drgState=1:3


function [sp_Orth,sp_Ret,TimeVars]=getFR_indcells(whichDatS,setParam)
% Time window for samples
Twin   = setParam.Twin;

%length of time to sample before stim
SponSampLength = setParam.SponSampLength;

WinOvrLap=setParam.Toverlap; %0=disjoint windows, 1=overlapping windows

if(Twin==2 || WinOvrLap==0) %if Twin=2, no choice, largest length of evoked window
    numEvok = (2/Twin); %assuming non-overlapping windows
else
    numEvok = 2*(2/Twin)-1; 
end

% Shift each trial by 15 seconds for display and calculation
StimShift = 15;

% To take samples of evoked activity
FirstEvok = round(StimShift/Twin) + 1;
LastEvok  = FirstEvok+numEvok-1;

% To take samples of spontaneous activity
LastSpon  = FirstEvok-1;
FirstSpon = round((StimShift-SponSampLength)/Twin)+1 ;   %SponSampLength (5) seconds before

TimeVars=[];
TimeVars.numEvok=numEvok;
TimeVars.StimShift=StimShift;
TimeVars.FirstEvok=FirstEvok;
TimeVars.LastEvok=LastEvok;
TimeVars.FirstSpon=FirstSpon;
TimeVars.LastSpon=LastSpon;

% Do all the files in the list!
nFiles = length(whichDatS);
dat_to_process = 1:nFiles;

%Output
sp_Ret=[];
sp_Orth=[];

for j1=dat_to_process
    % run processDataSelectShift.m
    [stim_trials,nEpoch,nOB,~]=processDataSelectShift(Twin,0,StimShift,whichDatS{j1}); %always set WinOvrLap=0 in processDataSelectShift b/c extractCountsOvrLap.m is wrong
    %stim_trials{1} is null? should not be used?
%     lenTrials = arrayfun(@(x) length(stim_trials{x}.OBcounts), 1:nEpoch-1);
    lenTrial = length(0:Twin:30.9);
    sp_RetTemp = zeros(lenTrial,floor(nEpoch/2),nOB);
    sp_OrthTemp = zeros(lenTrial,floor(nEpoch/2),nOB);
    for k=1:nOB
        OBcount_allTrial = zeros(lenTrial,nEpoch);
        how_many_stim_trials = 0;
        for k1=1:nEpoch
            if (stim_trials{k1}.stim_at_startT==1)
               % Is there a stimulus during this epoch??
               how_many_stim_trials = how_many_stim_trials+1;
               
               %Neuron ID value
%                if how_many_stim_trials==1
%                    nID=unique(stim_trials{k1}.OBuid);
%                    if(max(nID)~=nOB)
%                         keyboard
%                        nOB=min([max(nID) nOB]);
%                    end
%                end

                % Read data, put into array
                [nwin,nU]=size(stim_trials{k1}.OBcounts);

                if (k<=nU)   % Otherwise, No spikes belonging to this unit on that trial
                             % Recorded array was set based on max UID that
                             % appeared in that trial
                    if(nwin<lenTrial)
                        OBcount_allTrial(1:nwin,how_many_stim_trials)=stim_trials{k1}.OBcounts(1:nwin,k);
                    else
                        OBcount_allTrial(1:lenTrial,how_many_stim_trials)=stim_trials{k1}.OBcounts(1:lenTrial,k);
                    end
                end
            end
        end
        sp_RetTemp(:,:,k)=OBcount_allTrial(:,11:20);
        sp_OrthTemp(:,:,k)=OBcount_allTrial(:,1:10);
    end
    sp_Ret=cat(3,sp_Ret,sp_RetTemp);
    sp_Orth=cat(3,sp_Orth,sp_OrthTemp);
end
end