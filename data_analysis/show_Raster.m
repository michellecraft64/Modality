
% script to get a Raster, from no Drug
% Parameters
odorName='EB';

ind=[1 2 6 7 8 9 10 11]; %actual good rats
drugName='NoDrug';

%for Plots
tme=(-2:0.1:2); %must match with lines 17,19; length is =LastEvok-FirstSpon+1
onz=ones(1,length(tme));
for i=ind
        fileName=sprintf('Rat%d_IndCell_%s_%s.mat',i,odorName,drugName);
        %Load file, set parameters
        load(fileName) %Ortho/Retro raw spike counts (lenTime, numTrials, nID)
        numTrials=size(sOR,2)+size(sRET,2); %per or/ret cell
        numEvok=TimeVars.numEvok;
        StimShift=TimeVars.StimShift;
        FirstEvok=TimeVars.FirstEvok;
        LastEvok=FirstEvok+19; %TimeVars.LastEvok; %HARD CODED for opt twin-NEEDS CHANGING IF TWIN CHANGED
        LastSpon=TimeVars.LastSpon;
        FirstSpon=LastSpon-20; %TimeVars.FirstSpon; %HARD CODED for opt twin-NEEDS CHANGING IF TWIN CHANGED
        lenEvok=length(FirstEvok:1:LastEvok);
        lenSpon=length(FirstSpon:1:LastSpon);
        nOB=size(sOR,3);
        sRetTmp=sRET(FirstSpon:LastEvok,:,:);
        sOrTmp=sOR(FirstSpon:LastEvok,:,:);
        
        for clid=1:size(sOrTmp,3)
            figure
            hold on
            for trls=1:10
                plot(tme(sOrTmp(:,trls,clid)'==1),trls*onz(sOrTmp(:,trls,clid)'==1),'b.','MarkerSize',18)
            end
            for trls=11:20
                plot(tme(sRetTmp(:,trls-10,clid)'==1),trls*onz(sRetTmp(:,trls-10,clid)'==1),'r.','MarkerSize',18)
            end
            set(gca, 'Ydir', 'reverse')
            pause
        end
        close all
    end %all individual rats