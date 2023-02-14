function [stim_trials,nEpoch,nOB,nDat]=processDataSelectShift(Twin,WinOvrLap,StimShift,whichDatS)
%used in get_firing_rates.m and 
%relies on exclude_bad_units.m, extract_counts.m

% INPUTS
%   StimShift : Start each trial %% seconds BEFORE stimtime
%   whichDatS : EITHER an integer flag (OLD WAY)
%               or a file name (BETTER!)
whichDatS

if (ischar(whichDatS))
    infile = whichDatS;
    dSumm=[];
else
    if(whichDatS==0)
        
        infile = '../experimental_data/data012516_exp1+2';
        dSumm='../experimental_data/rat012506_1+2_summary.xlsx';
        
    elseif (whichDatS==1)
        infile = '../experimental_data/data012516_exp4+3';
        dSumm='../experimental_data/rat012506_4+3_summary.xlsx';
    elseif (whichDatS==3)
        infile = '../experimental_data/data040515_exp1+2';
        dSumm=[];
    elseif (whichDatS==4)
        infile = '../experimental_data/data040515_exp3+4';
        dSumm=[];
    end
end

%load summary of data
if (~isempty(dSumm))
    [nDat,tField]=xlsread(dSumm);  %only using cols 4(orth score), 5 (retro score), 7 (ortho p), 8 (retro p)
else
    nDat = [];
end

load(infile); %load spike stats 

%% Check file format here. If new format (Dec. 2019), translate.
if (~exist('OBspiketimes','var'))
    ind1 = find(spikes(:,3)==1);
    OBspiketimes = spikes(ind1,1);
    OBuid = uint32(spikes(ind1,2));
    stimtimes = expinfo(:,1)';
end

maxT  = ceil(max(OBspiketimes)); %no aPC

nOB  = max(OBuid);


%%%%% Firing rates
OB_fr  = zeros(nOB,1);

for k=1:nOB
    OB_fr(k)  = sum(OBuid==k)/maxT;
end

avgfr_OB = sum(OB_fr)/double(nOB);

% Check:
%disp(['Average firing rate of OB: ' sum(OB_fr)/double(nOB)]) 

% Decide on thresholds for exclusion
valid_fr_LB  = 5./(maxT-stimtimes(1));
valid_fr_UB  = 49.;

% dont count anything before stimtimes(1); throwout rates\notin[valid_fr_LB,valid_fr_UB]
%i_aftFirstStimOB=OBspiketimes>stimtimes(1); %logical vector
[OBspiketimes, OBuid, OB_unit_map]=exclude_bad_units(OBspiketimes,OBuid,...
    valid_fr_LB,valid_fr_UB,maxT);
nOB = length(OB_unit_map);

if (~isempty(nDat))
    %____ map the Summary data here, remove rows of units excluded ___
    nDat=nDat(OB_unit_map(:,1),:);
end

%stimtimes
%%%% Chop up data into trials
time_epochs  = [stimtimes-StimShift maxT];


if (time_epochs(1) > 0)
    % Add on "0", for period before shift
    time_epochs = [0 time_epochs];
% else
%     % Reset first time to zero
%     time_epochs(1) = 0;
end

%time_epochs

epoch_length = diff(time_epochs); 
nEpoch =length(epoch_length);

% Will include period before first stim
stim_trials = cell(nEpoch,1);
%PC_trials = cell(nEpoch,1);

%% Order the spikes
[OBsp_ord,I_OB]=sort(OBspiketimes,'ascend');
OBuid_ord = OBuid(I_OB);

%% First set up a structure
for k=1:nEpoch
    stim_trials{k} = struct('OBspikes',[],'OBuid',[],...
                            'startT',0,'endT',0,'stim_at_startT',0,'stim_Time_after_startT',0);
end

%% Enter some basic information
for k=1:nEpoch
    stim_trials{k}.startT = time_epochs(k);
    stim_trials{k}.endT   = time_epochs(k+1);
    
    jind = find(stimtimes >= stim_trials{k}.startT & stimtimes < stim_trials{k}.endT);
    if (~isempty(jind))
        stim_trials{k}.stim_at_startT = 1;    % Already set to zero by default
        
        stim_trials{k}.stim_Time_after_startT = stimtimes(jind)-stim_trials{k}.startT;
    end
end
        

%% OB spike times
lastind = 0;
for k=1:nEpoch
   firstind = lastind+1; 
   if (k==nEpoch)
       lastind = length(OBsp_ord);
   else
       lastind  = find(OBsp_ord > time_epochs(k+1),1)-1;
   end
   %[firstind lastind]
   
   % Record spikes relative to epoch start
   stim_trials{k}.OBspikes = OBsp_ord(firstind:lastind)-stim_trials{k}.startT;
   stim_trials{k}.OBuid    = OBuid_ord(firstind:lastind);
end

%% Spike counts
if(WinOvrLap==1)
    for k=1:nEpoch
        [nsp,cnts,Xcnts] = extractCountsOvrlap(stim_trials{k}.OBspikes,stim_trials{k}.OBuid,Twin);
       stim_trials{k}.OBcounts = cnts;
    end
else %original; disjoint windows
    for k=1:nEpoch
        [nsp,cnts,Xcnts] = extract_counts(stim_trials{k}.OBspikes,stim_trials{k}.OBuid,Twin);
        % Check counts
        cksum = sum(abs(nsp-(sum(cnts)' + Xcnts')));
        if (~cksum)
            % Ok!
            stim_trials{k}.OBcounts = cnts;
        end
    end
end



