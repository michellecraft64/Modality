function [A_M,A_P,synE,synI,w_MP,Atrn_M,Atrn_P,w_MPt]=...
            WC_EIcells(CpldWx,TimeVars,TimeConstants,xMin,InputCurrent,crr)
% Wilson-Cowen type rate model of network of 2 cells (MC - E, PGC - I)
% # glom = # MC, one-to-one #MC to PGC, there are NO GCs. 
% add alph fcn syn, 2 ODEs for each synE, synI
% approximated ODEs using RK4
% including plasticity on E->I; BUT changing A_p to A_m (pre-syn dependent)
% w_MP' = 1/tw*[ -w_MP + F(A_m) ] where F(x)=2./(1+exp((x-0)/1))-1

    Fp=@(x)(1./(1+exp((x-8)./1))-1);
    gnIsy=1;  %gain for I->E syn
    
    %Time Vars 
    dt=TimeVars.dt;
    lenSpon=TimeVars.lenSpon;
    lenEvok=TimeVars.lenEvok; %opt Twin = 900ms evoked (2s spont)
    lenTrnS=round(TimeVars.lenTrn/dt); %transient time
    tme=(0:dt:lenSpon+lenEvok)'; %time after transients from 0-tEnd; !!!ONLY used to get lenTime (dont need to start at -lenSpon
    lenTime = length(tme); %ms
    % Coupling strengths
    w_I=CpldWx.w_I; w_P=CpldWx.w_P;
    w_MPt=zeros(lenTrnS,1); wmp0=CpldWx.wmp0;
    % Time constants
    tau_M=TimeConstants.tau_M;
    tau_P=TimeConstants.tau_P;
    tdE=TimeConstants.tdE;  trE=TimeConstants.trE;
    tdI=TimeConstants.tdI;  trI=TimeConstants.trI;
    tw=TimeConstants.tw;
    % Minimal current for firing
    xM_min=xMin.xM;
    xP_min=xMin.xP;
    % Input Current
    I_M=InputCurrent.Mc; Itrn_M=I_M(1);
    I_P=InputCurrent.Pg; Itrn_P=I_P(1);
    % Initial activity for transients
    Atrn_M=zeros(lenTrnS,1);
    Atrn_P=zeros(lenTrnS,1);
    synEtrn=zeros(lenTrnS,1);
    asye=0; %not tracking
    synItrn=zeros(lenTrnS,1);
    asyi=0; %not tracking

    % Iterate through transients
    for i = 2:lenTrnS
        %First right hand side (RHS) of RK4 AKA k1
        rhs1_Mt = dt/tau_M*(-Atrn_M(i-1) + (w_I*Itrn_M - w_P.*synItrn(i-1) > xM_min).*(w_I*Itrn_M - w_P.*synItrn(i-1)));
        rhs1_Pt = dt/tau_P*(-Atrn_P(i-1) + (w_I*Itrn_P + (wmp0+w_MPt(i-1)).*synEtrn(i-1) > xP_min).*(w_I*Itrn_P + (wmp0+w_MPt(i-1)).*synEtrn(i-1)));
        asye1 = dt*(-asye/trE+Atrn_M(i-1));
        synE1 = dt/tdE*(-synEtrn(i-1)+asye);
        asyi1 = dt*(-asyi/trI+gnIsy*Atrn_P(i-1));
        synI1 = dt/tdI*(-synItrn(i-1)+asyi);
        w_MPt1 = dt/tw*crr*(Fp(Atrn_M(i-1))-w_MPt(i-1)); %Depr presynaptic = MC postsynaptic = PGC

        %k2
        rhs2_Mt = dt/tau_M*(-(Atrn_M(i-1)+0.5*rhs1_Mt) + (w_I*Itrn_M - w_P.*(synItrn(i-1)+0.5*synI1) > xM_min).*(w_I*Itrn_M - w_P.*(synItrn(i-1)+0.5*synI1)));
        rhs2_Pt = dt/tau_P*(-(Atrn_P(i-1)+0.5*rhs1_Pt) + (w_I*Itrn_P + (wmp0+w_MPt(i-1)+0.5*w_MPt1).*(synEtrn(i-1)+0.5*synE1) > xP_min).*(w_I*Itrn_P + (wmp0+w_MPt(i-1)+0.5*w_MPt1).*(synEtrn(i-1)+0.5*synE1)));
        asye2 = dt*(-(asye+.5*asye1)/trE+Atrn_M(i-1)+.5*rhs1_Mt);
        synE2 = dt/tdE*(-(synEtrn(i-1)+.5*synE1)+asye+.5*asye1);
        asyi2 = dt*(-(asyi+.5*asyi1)/trI+gnIsy*(Atrn_P(i-1)+.5*rhs1_Pt));
        synI2 = dt/tdI*(-(synItrn(i-1)+.5*synI1)+asyi+.5*asyi1);
        w_MPt2 = dt/tw*crr*(Fp(Atrn_M(i-1)+.5*rhs1_Pt)-(w_MPt(i-1)+.5*w_MPt1));
            
        %k3
        rhs3_Mt = dt/tau_M*(-(Atrn_M(i-1)+0.5*rhs2_Mt) + (w_I*Itrn_M - w_P.*(synItrn(i-1)+0.5*synI2) > xM_min).*(w_I*Itrn_M - w_P.*(synItrn(i-1)+0.5*synI2)));
        rhs3_Pt = dt/tau_P*(-(Atrn_P(i-1)+0.5*rhs2_Pt) + (w_I*Itrn_P + (wmp0+w_MPt(i-1)+0.5*w_MPt2).*(synEtrn(i-1)+0.5*synE2) > xP_min).*(w_I*Itrn_P + (wmp0+w_MPt(i-1)+0.5*w_MPt2).*(synEtrn(i-1)+0.5*synE2)));
        asye3 = dt*(-(asye+.5*asye2)/trE+Atrn_M(i-1)+.5*rhs2_Mt);
        synE3 = dt/tdE*(-(synEtrn(i-1)+.5*synE2)+asye+.5*asye2);
        asyi3 = dt*(-(asyi+.5*asyi2)/trI+gnIsy*(Atrn_P(i-1)+.5*rhs2_Pt));
        synI3 = dt/tdI*(-(synItrn(i-1)+.5*synI2)+asyi+.5*asyi2);
        w_MPt3 = dt/tw*crr*(Fp(Atrn_M(i-1)+.5*rhs2_Pt)-(w_MPt(i-1)+.5*w_MPt2));
        
        %k4
        rhs4_Mt = dt/tau_M*(-(Atrn_M(i-1)+rhs3_Mt) + (w_I*Itrn_M - w_P.*(synItrn(i-1)+synI3) > xM_min).*(w_I*Itrn_M - w_P.*(synItrn(i-1)+synI3)));
        rhs4_Pt = dt/tau_P*(-(Atrn_P(i-1)+rhs3_Pt) + (w_I*Itrn_P + (wmp0+w_MPt(i-1)+w_MPt3).*(synEtrn(i-1)+synE3) > xP_min).*(w_I*Itrn_P + (wmp0+w_MPt(i-1)+w_MPt3).*(synEtrn(i-1)+synE3)));
        asye4 = dt*(-(asye+asye3)/trE+Atrn_M(i-1)+rhs3_Mt);
        synE4 = dt/tdE*(-(synEtrn(i-1)+synE3)+asye+asye3);
        asyi4 = dt*(-(asyi+asyi3)/trI+gnIsy*(Atrn_P(i-1)+rhs3_Pt));
        synI4 = dt/tdI*(-(synItrn(i-1)+synI3)+asyi+asyi3);
        w_MPt4 = dt/tw*crr*(Fp(Atrn_M(i-1)+rhs3_Pt)-(w_MPt(i-1)+w_MPt3));
            
        %Final appx
        Atrn_M(i)=Atrn_M(i-1)+(1/6)*(rhs1_Mt+2*rhs2_Mt+2*rhs3_Mt+rhs4_Mt);
        Atrn_P(i)=Atrn_P(i-1)+(1/6)*(rhs1_Pt+2*rhs2_Pt+2*rhs3_Pt+rhs4_Pt);
        synEtrn(i)=synEtrn(i-1)+1/6*(synE1+2*synE2+2*synE3+synE4);
        synItrn(i)=synItrn(i-1)+1/6*(synI1+2*synI2+2*synI3+synI4);
        asye=asye+1/6*(asye1+2*asye2+2*asye3+asye4);
        asyi=asyi+1/6*(asyi1+2*asyi2+2*asyi3+asyi4);
        w_MPt(i)=w_MPt(i-1)+(1/6)*(w_MPt1+2*w_MPt2+2*w_MPt3+w_MPt4);
            
    end
    %Initialize output
    A_M=zeros(1,lenTime); A_M(1)=Atrn_M(end);
    A_P=zeros(1,lenTime); A_P(1)=Atrn_P(end);
    synE=zeros(1,lenTime); synE(1)=synEtrn(end);
    synI=zeros(1,lenTime); synI(1)=synItrn(end);
    w_MP=zeros(1,lenTime); w_MP(1)=w_MPt(end); 
    %keyboard
    %Sim ODE using RK4
    for j=2:lenTime
        %First right hand side (RHS) of RK4 AKA k1
        rhs1_Mt = dt/tau_M*(-A_M(j-1) + (w_I*I_M(j-1) - w_P.*synI(j-1) > xM_min).*(w_I*I_M(j-1) - w_P.*synI(j-1)));
        rhs1_Pt = dt/tau_P*(-A_P(j-1) + (w_I*I_P(j-1) + (wmp0+w_MP(j-1)).*synE(j-1) > xP_min).*(w_I*I_P(j-1) + (wmp0+w_MP(j-1)).*synE(j-1)));
        asye1 = dt*(-asye/trE+A_M(j-1));
        synE1 = dt/tdE*(-synE(j-1)+asye);
        asyi1 = dt*(-asyi/trI+gnIsy*A_P(j-1));
        synI1 = dt/tdI*(-synI(j-1)+asyi);
        w_MP1 = dt/tw*crr*(Fp(A_M(j-1))-w_MP(j-1)); %Depr presynaptic = MC postsynaptic = PGC
            
        %k2
        rhs2_Mt = dt/tau_M*(-(A_M(j-1)+0.5*rhs1_Mt) + (w_I*I_M(j-1) - w_P.*(synI(j-1)+0.5*synI1) > xM_min).*(w_I*I_M(j-1) - w_P.*(synI(j-1)+0.5*synI1)));
        rhs2_Pt = dt/tau_P*(-(A_P(j-1)+0.5*rhs1_Pt) + (w_I*I_P(j-1) + (wmp0+w_MP(j-1)+0.5*w_MP1).*(synE(j-1)+0.5*synE1) > xP_min).*(w_I*I_P(j-1) + (wmp0+w_MP(j-1)+0.5*w_MP1).*(synE(j-1)+0.5*synE1)));
        asye2 = dt*(-(asye+.5*asye1)/trE+A_M(j-1)+.5*rhs1_Mt);
        synE2 = dt/tdE*(-(synE(j-1)+.5*synE1)+asye+.5*asye1);
        asyi2 = dt*(-(asyi+.5*asyi1)/trI+gnIsy*(A_P(j-1)+.5*rhs1_Pt));
        synI2 = dt/tdI*(-(synI(j-1)+.5*synI1)+asyi+.5*asyi1);
        w_MP2 = dt/tw*crr*(Fp(A_M(j-1)+.5*rhs1_Pt)-(w_MP(j-1)+.5*w_MP1));
            
        %k3
        rhs3_Mt = dt/tau_M*(-(A_M(j-1)+0.5*rhs2_Mt) + (w_I*I_M(j-1) - w_P.*(synI(j-1)+0.5*synI2) > xM_min).*(w_I*I_M(j-1) - w_P.*(synI(j-1)+0.5*synI2)));
        rhs3_Pt = dt/tau_P*(-(A_P(j-1)+0.5*rhs2_Pt) + (w_I*I_P(j-1) + (wmp0+w_MP(j-1)+0.5*w_MP2).*(synE(j-1)+0.5*synE2) > xP_min).*(w_I*I_P(j-1) + (wmp0+w_MP(j-1)+0.5*w_MP2).*(synE(j-1)+0.5*synE2)));
        asye3 = dt*(-(asye+.5*asye2)/trE+A_M(j-1)+.5*rhs2_Mt);
        synE3 = dt/tdE*(-(synE(j-1)+.5*synE2)+asye+.5*asye2);
        asyi3 = dt*(-(asyi+.5*asyi2)/trI+gnIsy*(A_P(j-1)+.5*rhs2_Pt));
        synI3 = dt/tdI*(-(synI(j-1)+.5*synI2)+asyi+.5*asyi2);
        w_MP3 = dt/tw*crr*(Fp(A_M(j-1)+.5*rhs2_Pt)-(w_MP(j-1)+.5*w_MP2));
            
        %k4
        rhs4_Mt = dt/tau_M*(-(A_M(j-1)+rhs3_Mt) + (w_I*I_M(j-1) - w_P.*(synI(j-1)+synI3) > xM_min).*(w_I*I_M(j-1) - w_P.*(synI(j-1)+synI3)));
        rhs4_Pt = dt/tau_P*(-(A_P(j-1)+rhs3_Pt) + (w_I*I_P(j-1) + (wmp0+w_MP(j-1)+w_MP3).*(synE(j-1)+synE3) > xP_min).*(w_I*I_P(j-1) + (wmp0+w_MP(j-1)+w_MP3).*(synE(j-1)+synE3)));
        asye4 = dt*(-(asye+asye3)/trE+A_M(j-1)+rhs3_Mt);
        synE4 = dt/tdE*(-(synE(j-1)+synE3)+asye+asye3);
        asyi4 = dt*(-(asyi+asyi3)/trI+gnIsy*(A_P(j-1)+rhs3_Pt));
        synI4 = dt/tdI*(-(synI(j-1)+synI3)+asyi+asyi3);
        w_MP4 = dt/tw*crr*(Fp(A_M(j-1)+rhs3_Pt)-(w_MP(j-1)+w_MP3));
            
        %Final appx
        A_M(j)=A_M(j-1)+(1/6)*(rhs1_Mt+2*rhs2_Mt+2*rhs3_Mt+rhs4_Mt);
        A_P(j)=A_P(j-1)+(1/6)*(rhs1_Pt+2*rhs2_Pt+2*rhs3_Pt+rhs4_Pt);
        synE(j)=synE(j-1)+1/6*(synE1+2*synE2+2*synE3+synE4);
        synI(j)=synI(j-1)+1/6*(synI1+2*synI2+2*synI3+synI4);
        asye=asye+1/6*(asye1+2*asye2+2*asye3+asye4);
        asyi=asyi+1/6*(asyi1+2*asyi2+2*asyi3+asyi4);
        w_MP(j)=w_MP(j-1)+(1/6)*(w_MP1+2*w_MP2+2*w_MP3+w_MP4);
        
    end
end