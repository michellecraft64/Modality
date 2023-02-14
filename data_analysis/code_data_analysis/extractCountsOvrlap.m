function [ nspikes,counts, extracounts ] = extractCountsOvrlap( spikes, UIDs, winT )
% extractCountsOvrlap: turn spikes, UIDs into windowed spike counts
% very similar to extract_counts, except using overlapping windows (n+(n-1) total)
%
%   INPUTS:   spikes:     array of spikes
%             UIDs:       corresponding unit IDs
%             winT:       time window over which to count spikes
%
%   OUTPUTS:  counts:     array of counts
%             nspikes:    total spikes for each unit (can help to exclude
%                             non-spiking units)
%             extracounts:  list of spikes that occur after last winT
%             

%%% ASSUME minT = 0
maxT   = max(spikes);

maxID  = max(UIDs);


nwin  = 2*round(maxT/winT)-1; %n+(n-1) windows


counts   = zeros(nwin,maxID);
nspikes  = zeros(maxID,1);

% Keep counts after last time window
extracounts = zeros(1,maxID);



for k=1:maxID
    rightID    = find(UIDs==k);
    
    rightSpk   = spikes(rightID);
    if(length(rightSpk)>1) %only count spikes separated by 0.1 ms
        rightSpk = rightSpk([1==1;diff(rightSpk)>.0001]); 
    end
    
    %everything is counted twice besides 1st half bin & spikes after last half bins
    ind_cntTwice=rightSpk>=0.5*winT & rightSpk<maxT-.5*winT; %logical of spikes counted twice
    nspikes(k)=2*sum(ind_cntTwice)+sum(~ind_cntTwice); 
    
    if (~isempty(rightSpk))
        for k1=1:nwin
            %startT   = (k1-1)/2*winT; %never used b/c spikes are dynamically removed
            endT     = (k1+1)/2*winT;
            
            % First spike after end of interval
            lastind    = find(rightSpk >= endT-.5*winT,1); 
            
            %
            if (isempty(lastind))
                % Any remaining spikes are in THIS interval
                % Count them and move on to 
                counts(k1,k) = sum(rightSpk<endT);
                break;
            else
                % Any intervening spikes belong in this window
                counts(k1,k) = sum(rightSpk<endT);
                rightSpk = rightSpk(lastind:end); %remove spikes in [startT,startT+.5*winT)
            end
        end
        % Count any spikes AFTER last time window
        if (~isempty(rightSpk))
            extracounts(k) = length(find(rightSpk >= (nwin+1)/2*winT)); %no sliding windows at the end; just count
        end
         
    end
        
end

