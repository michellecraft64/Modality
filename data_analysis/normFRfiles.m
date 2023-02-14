% Script to create .mat files of normalized (spontaneous mean subtracted)
% individual cell FR for all drug preps
% Output: cell structure size 3 (number drug preps) x 2 (ortho/retro)
clear
% Parameters
odor = 2; 
odorName='EB';

% Initialize structure
NormStruct=cell(3,2);
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
    else
        disp('Param (drug_to_keep) needs to be set to 0, 1, or 2')
        return
    end
    numRats=length(ind);
    Zmn_Ret=[];
    Zmn_Orth=[];
    for i=ind
        fileName=sprintf('Rat%d_IndCell_%s_%s.mat',i,odorName,drugName);
        %Load file, set parameters
        load(fileName) %Ortho/Retro raw spike counts (lenTime, numTrials, nID)
        numTrials=size(sOR,2)+size(sRET,2); %per or/ret cell
        numEvok=TimeVars.numEvok;
        StimShift=TimeVars.StimShift;
        FirstEvok=TimeVars.FirstEvok;
        LastEvok=FirstEvok+8; %TimeVars.LastEvok; %HARD CODED for opt twin-NEEDS CHANGING IF TWIN CHANGED
        LastSpon=TimeVars.LastSpon;
        FirstSpon=LastSpon-19; %TimeVars.FirstSpon; %HARD CODED for opt twin-NEEDS CHANGING IF TWIN CHANGED
        lenEvok=length(FirstEvok:1:LastEvok);
        lenSpon=length(FirstSpon:1:LastSpon);
        nOB=size(sOR,3);
        if nOB~=size(sRET,3)
            keyboard
        end
        for k=1:nOB
            %Calculate evoked spike scores
            %Spontaneous mean subtract retro
            tmpEvkRet = sum(sRET(FirstEvok:LastEvok,:,k))/lenEvok;
            tmpSpnRet = sRET(FirstSpon:LastSpon,:,k);
            tmpZmn_Ret = tmpEvkRet-sum(tmpSpnRet)/lenSpon;
            %Spontaneous mean subtract ortho
            tmpEvkOrth = sum(sOR(FirstEvok:LastEvok,:,k))/lenEvok;
            tmpSpnOrth = sOR(FirstSpon:LastSpon,:,k);
            tmpZmn_Orth = tmpEvkOrth-sum(tmpSpnOrth)/lenSpon;
            %Append arrays
            Zmn_Orth=[Zmn_Orth;tmpZmn_Orth];
            Zmn_Ret=[Zmn_Ret;tmpZmn_Ret];
        end %all individual cells
    end %all individual rats
    NormStruct{drug_to_keep+1,1}=Zmn_Orth;
    NormStruct{drug_to_keep+1,2}=Zmn_Ret;
end %all drugs
fileName=sprintf('IndCell_normFR_%s.mat',odorName);
save(fileName,'NormStruct')