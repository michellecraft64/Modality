%create just total Decoding Accuracy mat file: DecodAcc.mat

load IndCellThrsh_OptWin_EB_2sSpon.mat

DecA_ND=[]; DecA_Bic=[]; DecA_Mus=[];

indWind=9; %index from 1,..,10; using 9=900ms evoked b/c largest p-val

for j=1:8 %total # rats, for ND
    DecA_ND=[DecA_ND; optThrsh_ND{j,4}(indWind,:).'];
    if j<=3 % bic & mus
        DecA_Bic=[DecA_Bic; optThrsh_Bic{j,4}(indWind,:).'];
        DecA_Mus=[DecA_Mus; optThrsh_Mus{j,4}(indWind,:).'];
    end
    if j==4 % bic
        DecA_Bic=[DecA_Bic; optThrsh_Bic{j,4}(indWind,:).'];
    end
end

save DecodAcc DecA_ND DecA_Bic DecA_Mus