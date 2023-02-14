function [metaData] = getMetaData_OdorDrug_SngleRat_fn(datfilestr,prestr)
% Get meta data for a single rat
%
% Now based on data in file: "expinfo" 
% Works for cleaned data: Dec. 2019
%
% 'expinfo' is a 20x3 matrix.  Each of the 20 rows corresponds to one olfactory 
% stimulus.  The first 10 stimuli are orthonasal.  The last 10 stimuli 
% are retronasal.  
% Column 1 contains the time of the stimulus onset in sec.  
% Column 2 contains the odor identity (1 = 1Hexanol, 2 = Ethyl Butyrate).  
% Column 3 indicates which (if any) pharmacological manipulation was done 
% to OB (0 = no drug, 1 = bicuculine, 2 = muscimol).

metaData = zeros(length(datfilestr),3);

for j1=1:length(datfilestr)
    load([prestr '/' datfilestr(j1).name],'-mat','expinfo');
    
    % 1st column won't be used anymore
    
    % 2nd column: odor ID
    metaData(j1,2) = expinfo(1,2);
    
    avgOdor = sum(expinfo(:,2))/20;
    if (round(avgOdor)~=avgOdor)
       warning(sprintf('Odor is not reported equal in %s', datfilestr.name));
    end
    
    %3rd column: drug/no drug
    metaData(j1,3) = expinfo(1,3);
    
    avgDrug = sum(expinfo(:,3))/20;
    if (round(avgDrug)~=avgDrug)
       warning(sprintf('Drug is not reported equal in %s', datfilestr.name));
    end
end