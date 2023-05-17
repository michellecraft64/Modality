% script to calc & save D (decode accur) for many values surrounding 
% fit to means (from dFitMeans.mat)
% dFitMeans was created from:
% >> save dFitMeans Am_O Am_R Ap_* InputCurrent_* Tmev alphMax* w_Parr dt sc aspPv Ntr idZer idEvD


save_flag=1; %if this is =1, WILL OVERWRITE EXISTING dJustDeAcc.mat

load dFitMeans.mat

%perturb from FitMeans
epsO=(-2: .25 : 4)'; 
epsR=(-6 : .5 : 9)';

lneO=length(epsO);
lneR=length(epsR);

DeAcc=zeros(lneO,lneR,3);

nuOint=sum(Am_O(:,idZer:idEvD),2)*dt/(1000*0.9);
nuRint=sum(Am_R(:,idZer:idEvD),2)*dt/(1000*0.9);

%initialize indices
indO=0; indR=0; k=0; prct=0;

for ind=1:(lneO*lneR*3) 

    [indO,indR,k]=ind2sub([lneO lneR 3],ind);
    %k: 1=ND, 2=Mus, 3=Bic
    
    meanOc=nuOint(k) + epsO(indO); %perturb from FitMeans
    rO=meanOc*sc(k,1);
    meanRc=nuRint(k) + epsR(indR); %perturb from FitMeans
    rR=meanRc*sc(k,2);

% Negative Binomial
sumSpks_O_Close=nbinrnd(rO,rO/(meanOc+rO),Ntr,1);
sumSpks_R_Close=nbinrnd(rR,rR/(meanRc+rR),Ntr,1);
%Set decode acc params
allRespn=[sumSpks_O_Close;sumSpks_R_Close];
mnOr=mean(sumSpks_O_Close,1);
mnRet=mean(sumSpks_R_Close,1);
if(mnRet>mnOr)
    lwThr=prctile(sumSpks_O_Close,3); %3rd percentile of smaller
    hiThr=prctile(sumSpks_R_Close,97); %97th percentile of larger
else
    lwThr=prctile(sumSpks_R_Close,3); %3rd percentile of smaller
    hiThr=prctile(sumSpks_O_Close,97); %97th percentile of larger
end
indT=(allRespn>=lwThr)&(allRespn<=hiThr);
thrsh=allRespn(indT); %consider subset for threshold values
Nthrs=length(thrsh);
RetAcc=zeros(Nthrs,1);
OrAcc=zeros(Nthrs,1);
ttlAcc=zeros(Nthrs,1);
for i=1:Nthrs %consider subset for threshold values
    Respn=ones(Ntr*2,1);
    if mnOr>mnRet
        Respn(allRespn>thrsh(i))=0;
    elseif mnRet>mnOr
        Respn(allRespn<thrsh(i))=0;
    else
        Respn(allRespn<thrsh(i))=0;
    end
    OrAcc(i)=1-sum(Respn(1:Ntr))/Ntr;
    RetAcc(i)=sum(Respn(Ntr+1:end))/Ntr;
    ttlAcc(i)=(OrAcc(i)+RetAcc(i))/2;
end % ---- end of 'new code' ----
[~,OptInd]=max(ttlAcc);

DeAcc(indO,indR,k)=ttlAcc(OptInd); %save results

if(save_flag==1 &&  indO==lneO) %save periodically
    save([pwd,'/dJustDeAcc'],'DeAcc','epsO','epsR')
    prct=ind/(lneO*lneR*3)
end

end


