% Script to do "off-response" analysis, MCs that increase firing after odor is removed 
% off-response: firing in (timOd_off,tim_end) is larger than firing in (0,timOd_off]

clear

timOd_off=1;  % time odor is removed, 1 sec after
tim_end=2;     % time ending to measure 'Off Response'

% Parameters
odorName='EB';

% Initialize structure
DrugStruct=cell(3,2); %Rows 1) Or 2) Ret / Columns 1) ND 2) Bic 3) Mus
% Loop over all drugs
for drug_to_keep = 0:2 % DRUG: no drug = 0; Bicu = 1 (less GABA_a inhib); Musc = 2 (more GABA_a inhib)
    %Load good rat data
    if drug_to_keep==0
        ind=[1 2 6 7 8 9 10 11]; %actual good rats
        drugName='NoDrug';
    elseif drug_to_keep==1
        ind=[1 2 8 11]; %actual good rats
        drugName='Bic';
    elseif drug_to_keep == 2
        ind=[6 9 10]; %actual good rats
        drugName='Mus';
    end
    numRats=length(ind);
    AllRat_sRET=[];
    AllRat_sOR=[];
    for i=ind
        fileName=sprintf('Rat%d_IndCell_%s_%s.mat',i,odorName,drugName);
        %Load file, set parameters
        load(fileName) %Ortho/Retro raw spike counts (lenTime, numTrials, nID)
        numTrials=size(sOR,2)+size(sRET,2); %per or/ret cell
        numEvok=TimeVars.numEvok;
        StimShift=TimeVars.StimShift;
        FirstEvok=TimeVars.FirstEvok;
        LastEvok=FirstEvok+49; %TimeVars.LastEvok; %HARD CODED for opt twin-NEEDS CHANGING IF TWIN CHANGED
        LastSpon=TimeVars.LastSpon;
        FirstSpon=LastSpon-20; %TimeVars.FirstSpon; %HARD CODED for opt twin-NEEDS CHANGING IF TWIN CHANGED
        lenEvok=length(FirstEvok:1:LastEvok);
        lenSpon=length(FirstSpon:1:LastSpon);
        nOB=size(sOR,3);
        sRetTmp=sRET(FirstSpon:LastEvok,:,:);
        sOrTmp=sOR(FirstSpon:LastEvok,:,:);
        AllRat_sRET=cat(3,AllRat_sRET,sRetTmp);
        AllRat_sOR=cat(3,AllRat_sOR,sOrTmp);
    end %all individual rats
    DrugStruct{drug_to_keep+1,1}=AllRat_sOR;
    DrugStruct{drug_to_keep+1,2}=AllRat_sRET;
end %all drugs
PSTH_Orth_ND=DrugStruct{1,1};  %counts, NOT PSTH (divide by time below)
PSTH_Ret_ND=DrugStruct{1,2};
PSTH_Orth_Bic=DrugStruct{2,1};
PSTH_Ret_Bic=DrugStruct{2,2};
PSTH_Orth_Mus=DrugStruct{3,1};
PSTH_Ret_Mus=DrugStruct{3,2};
% calc PSTH, avg over trials
mnPSTH_Ret_ND=squeeze(mean(PSTH_Ret_ND,2)); mnPSTH_Or_ND=squeeze(mean(PSTH_Orth_ND,2));
mnPSTH_Ret_Bic=squeeze(mean(PSTH_Ret_Bic,2)); mnPSTH_Or_Bic=squeeze(mean(PSTH_Orth_Bic,2));
mnPSTH_Ret_Mus=squeeze(mean(PSTH_Ret_Mus,2)); mnPSTH_Or_Mus=squeeze(mean(PSTH_Orth_Mus,2));

%% processing after get data

numND=size(mnPSTH_Or_ND,2); numBic=size(mnPSTH_Or_Bic,2); numMus=size(mnPSTH_Or_Mus,2);

% set times
tme=(-2:setParam.Twin:5)'; %must match with lines 32,34; length is =LastEvok-FirstSpon+1
idSt=(0-tme(1))/setParam.Twin+1; %index of tme=0, odor onset
idOdOff=(timOd_off-tme(1))/setParam.Twin+1; %index when odor removed
id_end=(tim_end-tme(1))/setParam.Twin+1; %index of end time

% fr_ND_or=mnPSTH_Or_ND./setParam.Twin; %if want Hz
% fr_Bic_or=mnPSTH_Or_Bic./setParam.Twin;
% fr_Mus_or=mnPSTH_Or_Mus./setParam.Twin;
% fr_ND_ret=mnPSTH_Ret_ND./setParam.Twin;
% fr_Bic_ret=mnPSTH_Ret_Bic./setParam.Twin;
% fr_Mus_ret=mnPSTH_Ret_Mus./setParam.Twin;

spk_odor_Or_ND=mean(mnPSTH_Or_ND(idSt:idOdOff,:))';
spk_after_Or_ND=mean(mnPSTH_Or_ND(idOdOff+1:id_end,:))';
perc_OffResp_Or_ND=sum(spk_odor_Or_ND<spk_after_Or_ND)/numND;
    spk_odor_Ret_ND=mean(mnPSTH_Ret_ND(idSt:idOdOff,:))';
    spk_after_Ret_ND=mean(mnPSTH_Ret_ND(idOdOff+1:id_end,:))';
    perc_OffResp_Ret_ND=sum(spk_odor_Ret_ND<spk_after_Ret_ND)/numND;

spk_odor_Or_Bic=mean(mnPSTH_Or_Bic(idSt:idOdOff,:))';
spk_after_Or_Bic=mean(mnPSTH_Or_Bic(idOdOff+1:id_end,:))';
perc_OffResp_Or_Bic=sum(spk_odor_Or_Bic<spk_after_Or_Bic)/numBic;
    spk_odor_Ret_Bic=mean(mnPSTH_Ret_Bic(idSt:idOdOff,:))';
    spk_after_Ret_Bic=mean(mnPSTH_Ret_Bic(idOdOff+1:id_end,:))';
    perc_OffResp_Ret_Bic=sum(spk_odor_Ret_Bic<spk_after_Ret_Bic)/numBic;
    
spk_odor_Or_Mus=mean(mnPSTH_Or_Mus(idSt:idOdOff,:))';
spk_after_Or_Mus=mean(mnPSTH_Or_Mus(idOdOff+1:id_end,:))';
perc_OffResp_Or_Mus=sum(spk_odor_Or_Mus<spk_after_Or_Mus)/numMus;
    spk_odor_Ret_Mus=mean(mnPSTH_Ret_Mus(idSt:idOdOff,:))';
    spk_after_Ret_Mus=mean(mnPSTH_Ret_Mus(idOdOff+1:id_end,:))';
    perc_OffResp_Ret_Mus=sum(spk_odor_Ret_Mus<spk_after_Ret_Mus)/numMus;

Modality={'Ortho';'Retro'};     
ND_OffCells=[perc_OffResp_Or_ND; perc_OffResp_Ret_ND];
Bic_OffCells=[perc_OffResp_Or_Bic; perc_OffResp_Ret_Bic];
Mus_OffCells=[perc_OffResp_Or_Mus; perc_OffResp_Ret_Mus];

T_EffSiz_DA=table(Modality,ND_OffCells,Bic_OffCells,Mus_OffCells)
