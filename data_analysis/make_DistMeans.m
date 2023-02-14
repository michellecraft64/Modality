%script to create DistMeans.mat file
clear
% Select parameters
odor = 2; %1 = Hex (food odor), 2 = EB (nonfood odor)
odorName='EB';

%Load good rat data
indCellFile=sprintf('IndCell_normFR_%s.mat',odorName);
load(indCellFile)
sumSpk_Orth_ND=NormStruct{1,1};
sumSpk_Ret_ND=NormStruct{1,2};
sumSpk_Orth_Bic=NormStruct{2,1};
sumSpk_Ret_Bic=NormStruct{2,2};
sumSpk_Orth_Mus=NormStruct{3,1};
sumSpk_Ret_Mus=NormStruct{3,2};
%calc mean spread
mnSpread_ND=abs(mean(sumSpk_Orth_ND,2)-mean(sumSpk_Ret_ND,2));
mnSpread_Bic=abs(mean(sumSpk_Orth_Bic,2)-mean(sumSpk_Ret_Bic,2));
mnSpread_Mus=abs(mean(sumSpk_Orth_Mus,2)-mean(sumSpk_Ret_Mus,2));

save DistMeans mnSpread_* 