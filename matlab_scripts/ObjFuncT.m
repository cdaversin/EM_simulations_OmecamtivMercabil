function [minObj,yt1,yp1]=ObjFuncT(om_params,points)

switch nargin
    case 1
        points = 1;
    otherwise

end

%initial conditions for state variables
    
    % Electrophysiology (41 variables):
    % X0=[v nai nass ki kss cai cass cansr cajsr m hf hs j hsp jp mL hL hLp a ...
    %     iF iS ap iFp iSp d ff fs fcaf fcas jca nca ffp fcafp xrf xrs xs1 xs2 xk1...
    %     Jrelnp Jrelp CaMKt]';
    %+Mechanical 
    %X0(42:48)=[XS; XW; CaTrpn; TmB; Zetas; Zetaw; Cd];

%%

CL=1000;%pacing cycle length in ms
beats=10; % enough if I only modify Land parameters
mech.cell='endo';


options=[];%options for ode solver
tspan=[0 CL]; %[0 CL] variable dt

mech.isacs = 0; % 0 / 1: Isacs
mech.emcoupling = 1; % 0: ORd only; 1: ORd-Land interaction and feedback
    mech.mode = 'intact'; % 'skinned'/'intact'
    mech.lambda = 1; %1.1 max in Isacs; default Land: 1
    mech.dLambda = 0;
mech.calib = 1; %1: Margara calibration in ORd_dutta (approx to ORdmm)

model=@model_Tor_Land_OM; 
% N endo SS 
    X0 = [-89.0831180378458,12.0590970262148,12.0594189998681,143.661564172166,143.661522837075,7.28595493404270e-05,6.34904342041025e-05,1.53802427606100,1.53936005481898,0.000752189976974158,0.681021330125326,0.834865861495171,0.834900340648339,0.834742315524641,0.000153311232111416,0.537306964102465,0.297408268761484,0.000930720896426404,0.999631110120037,0.596153805396994,0.000474211033721301,0.999631107510222,0.658870731010089,6.07929048437019e-33,0.999999993191107,0.939984414627503,0.999999993191090,0.999899545002554,0.999978532989767,0.000447862732201756,0.000762504648734618,0.999999993191719,0.999999993191876,0.240799807585714,0.000170792294271546,-2.45639084191538e-25,0.0112323694607491,0.998132573935102,0.000832009126478113,0.000686086625522053,0.000337594734372013,1.17336380851157e-05,1.65492398556338e-23,0.000156961427588819,0.000236320366283090,0.00811841967054613,0.999370262041957,0,0,0];

if points==1
    conc=[0 0.1 0.3 1 3];
else
    multx=-2:0.2:2; conc1=10.^multx; %curve
    conc2=[0.1 3]; %points for TTP
    conc=[conc1,conc2];
end


% Method (Hill-EC50)
    kuw_sf=1+om_params(1)./(1+(om_params(2)./conc).^om_params(3));
    kwu_sf=1-om_params(4)./(1+(om_params(5)./conc).^om_params(6));
    kws_sf=1-om_params(7)./(1+(om_params(8)./conc).^om_params(9));


for j=1:length(conc)
    
    mech.factor1=kuw_sf(j);  %OM
    mech.factor2=kwu_sf(j);
    mech.factor3=kws_sf(j);

    for n=1:beats

        [time X]=ode15s(@(t,y)model(t,y,mech,1),tspan,X0);
        X0=X(size(X,1),:);

    end


    %calculate and name dependent variables for the final beat in the
    %simulation (i.e. currents and fluxes)
    for i=[1:size(X,1)]
        %IsJs=model(time(i),X(i,:),mech,0);
        tension=model(time(i),X(i,:),mech,2);

        Ttot(i)=tension(1);
        Ta(i)=tension(2);
        Tp(i)=tension(3);
    end

    tens(j).time=time;
    tens(j).Ta=Ta;
    [rt50(j),tp(j),maxTa(j),minTa(j)]=Txx(time,Ta,50);
    
    clear T*

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rel_pTa=(maxTa-24.702)/24.702;
Atp=(tp-171.01)./171.01;

if points==1
    yexp1=[0 0.1404 0.4741 1.0046 1.2441];
    yt1=rel_pTa;

    yexp2=[0.1,0.5];
    yp2=Atp([2,4]);
else
    %Cond1)
    yexp1=1.33./(1+(0.47./conc1).^1.41); 
    yt1=rel_pTa(1:length(conc1));  %output:multiple points
    yp1=Atp(1:length(conc1));      %output:multiple points
    %Cond2)
    yexp2=[0.1,0.5];
    yp2=Atp(end-1:end);  %2 points to compare with exp
end


%--------------------------------------------------------------------------
% Mean squared error to minimize

cond1=mean((yexp1-yt1).^2);
cond2=mean((yexp2-yp2).^2);

minObj=cond1+10*cond2; %with different weights
%minObj=cond1+cond2;
%--------------------------------------------------------------------------


